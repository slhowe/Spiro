#define BIT_OFFSET 103
#define V_S 5
#define FUNCTION_PERIOD_US 3333 //300 Hz

#define MIN_FLOW_RATE 0.1 //L\s 
#define WAIT_DURATION 4000000 //us

#define BLINKY 2
#define SOLTOGGLE 17

/* Predefine functions
 */
void setup();
float abs_pressure_conversion();
void convert_to_Pa_ASDX(int* value, int Pmin, int Pmax);
void convert_to_Pa_RBIP(int* value, unsigned long Pa_per_bit);
//void convert_to_Pa_HCLA(int* value);
void convert_to_Pa_SSC(int* value);
static inline int8_t sgn(float val);
int convert_Pa_to_L_min(float value, float sqr_const, float lin_const);
float convert_Pa_to_L_s(float value);
/* Timer runs in micros
 * This program will run for ~70 mins before overflow occurs
 */
unsigned long function_run_time = 0;
unsigned long current_time = FUNCTION_PERIOD_US; //Function to run first loop
unsigned long shutter_time = 0;

int buffer[15];
bool occluding = false;

long occlusion_duration[6] = {300000, 200000, 150000, 100000, 75000, 50000}; //us
int occlusion_index = 0;


/* Main loop:
 *  Read analogue data
 *  Convert data to pressure
 *  Convert differential pressure to flow
 *  Print results to serial
 */
void loop() {
  // initial setup of serial, constants etc.
  setup();
  //static float Pa_per_bit = abs_pressure_conversion();
  static float sqr_const = -0.00047541; // From calibration
  static float lin_const = 0.4269203; // From calibration
  int high_pressure_measurement;
  int low_pressure_measurement;
  float flow;


  while(1){
    //data reading
    if(current_time >= (function_run_time + FUNCTION_PERIOD_US)){
      function_run_time = micros();


      // read the input
      high_pressure_measurement = analogRead(A10); //High pressure sensor
      low_pressure_measurement = analogRead(A9); //High res sensor

      // Convert reading to pressure
      convert_to_Pa_SSC(&low_pressure_measurement);
      convert_to_Pa_ASDX(&high_pressure_measurement);
        
      /*
       * NOTE
       * Expiration is negative flow
       * Check connection is correct way around!
       */
      low_pressure_measurement *= -1;
      high_pressure_measurement -= 33; //adjust for error
      
      flow = convert_Pa_to_L_s(low_pressure_measurement);
      
      // print out the values
      Serial.print(low_pressure_measurement); //High res sensor
      Serial.print(",");
      Serial.print(flow);
      Serial.print(",");
      Serial.print(-high_pressure_measurement);
      Serial.print(",");
      Serial.print(function_run_time);
      Serial.print("\n");
    }

    //shutter re-opening
    if(occluding && current_time >= (shutter_time + occlusion_duration[occlusion_index])){
      shutter_time = micros();
      digitalWrite(SOLTOGGLE, 0);
      digitalWrite(BLINKY, !digitalRead(BLINKY));
      occluding = false;
      occlusion_index++;
      if (occlusion_index >= 6){
        occlusion_index = 0;
      }
    }
    // shutter closing timing
    else if(!occluding && current_time >= (shutter_time + WAIT_DURATION)){
      shutter_time = micros();
      digitalWrite(SOLTOGGLE, 1);
      digitalWrite(BLINKY, !digitalRead(BLINKY));
      occluding = true;
    }
    
    current_time = micros();    
    
  }
}

/* Set baud rate and open communications
 */
void setup() {
  // initialize serial communication at 115200 bits per second:
  Serial.begin(115200);  
  /*
   * A9 Low pressure sensor
   * A10 High pressure sensor
   * A8 solenoid voltage
   * 17 solenoid toggle
   * 2 blinky
   */
  pinMode(2, OUTPUT); //BLINKY
  pinMode(17, OUTPUT); //solToggle
  
}

/* The conversion constant from bits to Pa
 * for the absolute pressure sensor.
 */
float abs_pressure_conversion(void){
  // Values from datasheet (R B I P 0 0 1 D U)
  float Pa_per_mV = 1.72369;
  float mV_per_bit = 4.887585;
  float conversion = Pa_per_mV * mV_per_bit;
  return conversion;
}

/* Convert bit value returned by analogue pin
 * to Pa for the differential pressure sensor 
 * (A S D X R R X 0 0 5 P D A A 5).
 * Equation is (datasheet --> A type cal)
 *  (Vmeas - 0.1*Vs)(Pmax - Pmin) + Pmin = Pmeas
 *  -----------------------------
 *           (0.8*Vs)
 */
void convert_to_Pa_ASDX(int* value){
  /* Convert from bits to volts
   * 4.887585 mv per bit
   * 1000 mv per V
   */
  static float Pmax = 27579; //Pa
  static float Pmin = -27579; //Pa
  
  float temp_value = *value * 4.887585 / 1000.0;
  
  temp_value -= 0.1*V_S;
  temp_value *= (Pmax - Pmin);
  temp_value /= 0.8*V_S;
  temp_value += Pmin;
  
  *value = (int)temp_value;
}

/* Convert bit value returned by analogue pin
 * to Pa for the differential pressure sensor
 * (S S C S N B N 0 0 1 N D A A 5)
 */
void convert_to_Pa_SSC(int* value){
  static float Pmin = -248.84; //Pa
  static float Pmax = 248.84; //Pa
  
  // Convert bits to voltage
  float temp_value = *value * (5.0/1023.0);
  
  // linear equation from datasheet specifications
  temp_value -= 0.5; //0.1*Vs
  temp_value *= (Pmax - Pmin);
  temp_value /= 4; //0.8*Vs
  temp_value += Pmin;

  *value = (int)temp_value;
}

/* Get the sign of a float
 * -1 if negative
 *  0 if 0
 *  1 if positive
 */
static inline int8_t sgn(float val) {
 if (val < 0) return -1;
 if (val==0) return 0;
 return 1;
}

float convert_Pa_to_L_s(int value){
  float flow = 0;
  if (value >= 0) //inspiration
    //flow = 0.0071*value;i
    //adjust 50% for error in venturi test (Q^2)
    flow = 0.0088*value;
  else
    //flow = 0.0075*value;
    //adjust 50% for error in venturi test (Q^2)
    flow = 0.0092*value;
  return flow;
}

/* Convert pressure formn differential sensor
 * to flow. Uses emperic equation from calibration.
 */
int convert_Pa_to_L_min(int value, float sqr_const, float lin_const){
  int8_t sn = sgn(value);
  if (sn == 0) return 0;

  // old method
  float value_abs = abs(value);
  float quad = value_abs * value_abs * sqr_const;
  float lin = value_abs * lin_const;
  float flow = quad + lin; // Round to nearest int

  
  return flow;
}
