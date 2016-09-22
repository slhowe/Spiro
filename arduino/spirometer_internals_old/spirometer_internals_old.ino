#define mV_PER_Pa 8
#define FIVE_mV 5000
#define BIT_RANGE 1023
#define CONVERSION_OFFSET 1000
#define BIT_OFFSET 460 // 0 kPa in bit value

// function definitions here, have loop at top
void setup();
unsigned int pressure_conversion();
void convert_to_mV(unsigned int value);

void loop() {
  setup();
  static int Pa_per_bit = pressure_conversion();
  //static float sqr_const = -0.000475; // From calibration
  //static float lin_const = 0.427; // From calibration

  while(1){
    
    // read the input
    // Convert to voltage
    int diff_pressure = analogRead(A0);

    convert_to_Pa(&diff_pressure, Pa_per_bit);
    //int flow = convert_to_L_min((float)diff_pressure, sqr_const, lin_const);

    // print out the value you read
    Serial.println(diff_pressure);
    //Serial.print(",");
    //Serial.println(flow);

  }
}

void setup() {
  // initialize serial communication at 9600 bits per second:
  Serial.begin(9600);
}

unsigned int pressure_conversion(void){
  // convert from 1023 bits to Pa
  unsigned long conversion = FIVE_mV / mV_PER_Pa;
  conversion = conversion * CONVERSION_OFFSET;
  conversion = conversion / BIT_RANGE;
  return int(conversion);
}

void convert_to_Pa(int* value, unsigned long Pa_per_bit){
  long temp_value = *value - BIT_OFFSET;
  temp_value *= Pa_per_bit;
  temp_value /= CONVERSION_OFFSET;
  *value = int(temp_value);
}

static inline int8_t sgn(float val) {
 if (val < 0) return -1;
 if (val==0) return 0;
 return 1;
}

int convert_to_L_min(float value, float sqr_const, float lin_const){
  // Using Q = aP^2 + bP
  int8_t sn = sgn(value);
  if (sn == 0) return 0;
  
  value = abs(value);
  
  float quad = value * value * sqr_const;
  float lin = value * lin_const;
  float flow = quad + lin + 0.5; // Round to nearest int
  
  return (int)flow * sn * -1;
}

