#define BIT_OFFSET 103
#define V_S 5

/* Predefine functions
 */
void setup();
float abs_pressure_conversion();
void convert_to_Pa_ASDX(int* value, int Pmin, int Pmax);
void convert_to_Pa_RBIP(int* value, unsigned long Pa_per_bit);
static inline int8_t sgn(float val);
int convert_Pa_to_L_min(float value, float sqr_const, float lin_const);

/* Main loop:
 *  Read analogue data
 *  Convert data to pressure
 *  Convert differential pressure to flow
 *  Print results to serial
 */
void loop() {
  // initial setup of serial, constants etc.
  setup();
  static float Pa_per_bit = abs_pressure_conversion();
  static float Pmax = 6205; //Pa
  static float Pmin = -6205; //Pa
  static float sqr_const = -0.000475; // From calibration
  static float lin_const = 0.427; // From calibration

  while(1){
    // read the input
    int flow_measurement = analogRead(A1);
    int pressure_measurement = analogRead(A0);

    // Convert reading to pressure
    convert_to_Pa_ASDX(&pressure_measurement, Pmin, Pmax);
    convert_to_Pa_RBIP(&flow_measurement, Pa_per_bit);
    flow_measurement += 4;

    // Convert differential pressure to flow
    int flow = convert_Pa_to_L_min((float)flow_measurement, sqr_const, lin_const);

    // Timestamp
    int timestamp = millis();

    // print out the values
    Serial.print(pressure_measurement);
    Serial.print(",");
    Serial.print(flow);
    Serial.print(",");
    Serial.print(timestamp);
    Serial.print("\n");

  }
}

/* Set baud rate and open communications
 */
void setup() {
  // initialize serial communication at 9600 bits per second:
  Serial.begin(9600);
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
 * (A S D X R R X 0 0 1 P D A A 5).
 * Equation is (datasheet --> A type cal)
 *  (Vmeas - 0.1*Vs)(Pmax - Pmin) + Pmin = Pmeas
 *  -----------------------------
 *           (0.8*Vs)
 */
void convert_to_Pa_ASDX(int* value, int Pmin, int Pmax){
  /* Convert from bits to volts
   * 4.887585 mv per bit
   * 1000 mv per V
   */
  float temp_value = *value * 4.887585 / 1000.0;
  
  temp_value -= 0.1*V_S;
  temp_value *= (Pmax - Pmin);
  temp_value /= 0.8*V_S;
  temp_value += Pmin;
  
  *value = (int)temp_value;
}

/* Convert bit value returned by analogue pin
 * to Pa for the absolute pressure sensor
 * (R B I P 0 0 1 D U)
 */
void convert_to_Pa_RBIP(int* value, unsigned long Pa_per_bit){
  // 0 Pa not at 0 bits
  float temp_value = *value - BIT_OFFSET;

  // convert bits to Pa
  temp_value *= Pa_per_bit;

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

/* Convert pressure formn differential sensor
 * to flow. Uses emperic equation from calibration.
 */
int convert_Pa_to_L_min(float value, float sqr_const, float lin_const){
  // Using Q = aP^2 + bP
  int8_t sn = sgn(value);
  if (sn == 0) return 0;
  
  value = abs(value);
  
  float quad = value * value * sqr_const;
  float lin = value * lin_const;
  float flow = quad + lin + 0.5; // Round to nearest int
  
  return (int)flow * sn * -1;
}

