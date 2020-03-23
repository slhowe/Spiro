EESchema Schematic File Version 4
EELAYER 30 0
EELAYER END
$Descr A4 11693 8268
encoding utf-8
Sheet 1 1
Title "Solenoid"
Date ""
Rev ""
Comp ""
Comment1 ""
Comment2 ""
Comment3 ""
Comment4 ""
$EndDescr
$Comp
L Device:CP1 Csuper1
U 1 1 5C50FC3F
P 5600 1200
F 0 "Csuper1" H 5715 1246 50  0000 L CNN
F 1 "1F" H 5715 1155 50  0000 L CNN
F 2 "SpiroParts:SCMR18H105PRBB0-SuperCap1F" H 5638 1050 50  0001 C CNN
F 3 "~" H 5600 1200 50  0001 C CNN
	1    5600 1200
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR0105
U 1 1 5C50FCBF
P 5600 1650
F 0 "#PWR0105" H 5600 1400 50  0001 C CNN
F 1 "GND" H 5605 1477 50  0000 C CNN
F 2 "" H 5600 1650 50  0001 C CNN
F 3 "" H 5600 1650 50  0001 C CNN
	1    5600 1650
	1    0    0    -1  
$EndComp
Text Label 5800 1050 0    50   ~ 0
solenoidPWR
Wire Wire Line
	5450 1050 5600 1050
Connection ~ 5600 1050
Wire Wire Line
	5600 1050 5800 1050
$Comp
L power:GND #PWR0107
U 1 1 5C52319F
P 2950 1400
F 0 "#PWR0107" H 2950 1150 50  0001 C CNN
F 1 "GND" H 2955 1227 50  0000 C CNN
F 2 "" H 2950 1400 50  0001 C CNN
F 3 "" H 2950 1400 50  0001 C CNN
	1    2950 1400
	1    0    0    -1  
$EndComp
Wire Wire Line
	2250 1000 2350 1000
Connection ~ 2350 1000
$Comp
L Device:C C0
U 1 1 5C5231A9
P 2350 1500
F 0 "C0" H 2465 1546 50  0000 L CNN
F 1 "1uF" H 2465 1455 50  0000 L CNN
F 2 "Capacitor_SMD:C_0603_1608Metric" H 2388 1350 50  0001 C CNN
F 3 "~" H 2350 1500 50  0001 C CNN
	1    2350 1500
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR0108
U 1 1 5C5231B8
P 3750 1650
F 0 "#PWR0108" H 3750 1400 50  0001 C CNN
F 1 "GND" H 3755 1477 50  0000 C CNN
F 2 "" H 3750 1650 50  0001 C CNN
F 3 "" H 3750 1650 50  0001 C CNN
	1    3750 1650
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR0109
U 1 1 5C5231BE
P 2350 1650
F 0 "#PWR0109" H 2350 1400 50  0001 C CNN
F 1 "GND" H 2355 1477 50  0000 C CNN
F 2 "" H 2350 1650 50  0001 C CNN
F 3 "" H 2350 1650 50  0001 C CNN
	1    2350 1650
	1    0    0    -1  
$EndComp
$Comp
L Device:C C2
U 1 1 5C5231C4
P 3750 1500
F 0 "C2" H 3865 1546 50  0000 L CNN
F 1 "2.2uF" H 3865 1455 50  0000 L CNN
F 2 "Capacitor_SMD:C_0603_1608Metric" H 3788 1350 50  0001 C CNN
F 3 "~" H 3750 1500 50  0001 C CNN
	1    3750 1500
	1    0    0    -1  
$EndComp
Text Label 2250 1000 2    50   ~ 0
Vin
Wire Notes Line
	1900 2000 1900 600 
Text Notes 3450 700  0    50   ~ 0
5V regulator
Wire Notes Line
	4600 600  4600 2000
$Comp
L Device:CP1 Csuper2
U 1 1 5C590AED
P 5600 1500
F 0 "Csuper2" H 5715 1546 50  0000 L CNN
F 1 "1F" H 5715 1455 50  0000 L CNN
F 2 "SpiroParts:SCMR18H105PRBB0-SuperCap1F" H 5638 1350 50  0001 C CNN
F 3 "~" H 5600 1500 50  0001 C CNN
	1    5600 1500
	1    0    0    -1  
$EndComp
Text Label 5450 1050 2    50   ~ 0
Vin
Text Label 6900 1350 2    50   ~ 0
solenoidToggle
$Comp
L Device:R R2
U 1 1 5C5193D8
P 7050 1350
F 0 "R2" V 6950 1350 50  0000 C CNN
F 1 "1K" V 7050 1350 50  0000 C CNN
F 2 "Resistor_SMD:R_0201_0603Metric" V 6980 1350 50  0001 C CNN
F 3 "~" H 7050 1350 50  0001 C CNN
	1    7050 1350
	0    1    1    0   
$EndComp
$Comp
L power:GND #PWR0106
U 1 1 5C5194EB
P 7500 1550
F 0 "#PWR0106" H 7500 1300 50  0001 C CNN
F 1 "GND" H 7505 1377 50  0000 C CNN
F 2 "" H 7500 1550 50  0001 C CNN
F 3 "" H 7500 1550 50  0001 C CNN
	1    7500 1550
	1    0    0    -1  
$EndComp
$Comp
L Connector:Screw_Terminal_01x02 J3
U 1 1 5C592068
P 6800 2050
F 0 "J3" H 6880 2042 50  0000 L CNN
F 1 "Screw_Terminal_01x02" H 6880 1951 50  0000 L CNN
F 2 "TerminalBlock_Phoenix:TerminalBlock_Phoenix_PT-1,5-2-5.0-H_1x02_P5.00mm_Horizontal" H 6800 2050 50  0001 C CNN
F 3 "~" H 6800 2050 50  0001 C CNN
	1    6800 2050
	1    0    0    -1  
$EndComp
Text Label 7500 1150 1    50   ~ 0
solenoidGND
Text Label 6600 2150 2    50   ~ 0
solenoidGND
Text Label 6600 2050 2    50   ~ 0
solenoidPWR
Text Label 5650 4150 0    50   ~ 0
PressureSensorLow1
Text Label 5650 4050 0    50   ~ 0
PressureSensorLow2
Text Label 5650 4250 0    50   ~ 0
PressureSensorHigh
Text Label 5650 3650 0    50   ~ 0
solenoidToggle
$Comp
L Connector:Barrel_Jack_Switch J1
U 1 1 5C5B5393
P 1000 1550
F 0 "J1" H 1055 1867 50  0000 C CNN
F 1 "Barrel_Jack_Switch" H 1055 1776 50  0000 C CNN
F 2 "SpiroParts:PJ-025-powerJack" H 1050 1510 50  0001 C CNN
F 3 "~" H 1050 1510 50  0001 C CNN
	1    1000 1550
	1    0    0    -1  
$EndComp
Text Label 1300 1450 0    50   ~ 0
Vin
$Comp
L power:GND #PWR0101
U 1 1 5C5B5512
P 1300 1650
F 0 "#PWR0101" H 1300 1400 50  0001 C CNN
F 1 "GND" H 1305 1477 50  0000 C CNN
F 2 "" H 1300 1650 50  0001 C CNN
F 3 "" H 1300 1650 50  0001 C CNN
	1    1300 1650
	1    0    0    -1  
$EndComp
$Comp
L solenoidParts:ASDXRRX001PDAA5 PSH1
U 1 1 5C5B7767
P 9250 950
F 0 "PSH1" H 9600 1000 50  0000 C CNN
F 1 "ASDXRRX001PDAA5" H 9450 -154 50  0000 C CNN
F 2 "SpiroParts:ASDXRRX001PDAA5-PressureSensor" H 9250 950 50  0001 C CNN
F 3 "" H 9250 950 50  0001 C CNN
	1    9250 950 
	1    0    0    -1  
$EndComp
$Comp
L solenoidParts:SSCSNBN001NDAA5 PSL1
U 1 1 5C5B788B
P 9200 2650
F 0 "PSL1" H 9600 2700 50  0000 C CNN
F 1 "SSCSNBN001NDAA5" H 9425 1946 50  0000 C CNN
F 2 "SpiroParts:SSCSNBN001NDAA5" H 9200 2650 50  0001 C CNN
F 3 "" H 9200 2650 50  0001 C CNN
	1    9200 2650
	1    0    0    -1  
$EndComp
$Comp
L solenoidParts:SSCSNBN001NDAA5 PSL2
U 1 1 5C5B7901
P 9200 3850
F 0 "PSL2" H 9600 3900 50  0000 C CNN
F 1 "SSCSNBN001NDAA5" H 9425 3146 50  0000 C CNN
F 2 "SpiroParts:SSCSNBN001NDAA5" H 9200 3850 50  0001 C CNN
F 3 "" H 9200 3850 50  0001 C CNN
	1    9200 3850
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR0102
U 1 1 5C5B8307
P 9500 1700
F 0 "#PWR0102" H 9500 1450 50  0001 C CNN
F 1 "GND" H 9505 1527 50  0000 C CNN
F 2 "" H 9500 1700 50  0001 C CNN
F 3 "" H 9500 1700 50  0001 C CNN
	1    9500 1700
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR0103
U 1 1 5C5B8350
P 9400 3000
F 0 "#PWR0103" H 9400 2750 50  0001 C CNN
F 1 "GND" H 9405 2827 50  0000 C CNN
F 2 "" H 9400 3000 50  0001 C CNN
F 3 "" H 9400 3000 50  0001 C CNN
	1    9400 3000
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR0104
U 1 1 5C5B8418
P 9400 4200
F 0 "#PWR0104" H 9400 3950 50  0001 C CNN
F 1 "GND" H 9405 4027 50  0000 C CNN
F 2 "" H 9400 4200 50  0001 C CNN
F 3 "" H 9400 4200 50  0001 C CNN
	1    9400 4200
	1    0    0    -1  
$EndComp
Text Label 9100 2750 2    50   ~ 0
VCC
Text Label 9100 3950 2    50   ~ 0
VCC
Text Label 9150 1000 2    50   ~ 0
VCC
Text Label 9750 1250 0    50   ~ 0
PressureSensorHigh
Text Label 9750 2750 0    50   ~ 0
PressureSensorLow1
Text Label 9750 3950 0    50   ~ 0
PressureSensorLow2
$Comp
L baseBoard-rescue:ATmega32U4-MU-MCU_Microchip_ATmega U4
U 1 1 5C5BC72D
P 5050 5150
F 0 "U4" H 5500 6900 50  0000 C CNN
F 1 "ATmega32U4-MU" H 5050 2850 50  0000 C CNN
F 2 "Package_DFN_QFN:QFN-44-1EP_7x7mm_P0.5mm_EP5.2x5.2mm" H 5050 5150 50  0001 C CIN
F 3 "http://ww1.microchip.com/downloads/en/DeviceDoc/Atmel-7766-8-bit-AVR-ATmega16U4-32U4_Datasheet.pdf" H 5050 5150 50  0001 C CNN
	1    5050 5150
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR0110
U 1 1 5C5BD6E7
P 5050 6950
F 0 "#PWR0110" H 5050 6700 50  0001 C CNN
F 1 "GND" H 5055 6777 50  0000 C CNN
F 2 "" H 5050 6950 50  0001 C CNN
F 3 "" H 5050 6950 50  0001 C CNN
	1    5050 6950
	1    0    0    -1  
$EndComp
Wire Wire Line
	4950 6950 5050 6950
Connection ~ 5050 6950
$Comp
L Connector:USB_B_Mini J2
U 1 1 5C5BED2E
P 1000 4100
F 0 "J2" H 1055 4567 50  0000 C CNN
F 1 "USB_B_Mini" H 1055 4476 50  0000 C CNN
F 2 "Connector_USB:USB_Mini-B_Lumberg_2486_01_Horizontal" H 1150 4050 50  0001 C CNN
F 3 "~" H 1150 4050 50  0001 C CNN
	1    1000 4100
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR0111
U 1 1 5C5BEDDE
P 1000 4500
F 0 "#PWR0111" H 1000 4250 50  0001 C CNN
F 1 "GND" H 1005 4327 50  0000 C CNN
F 2 "" H 1000 4500 50  0001 C CNN
F 3 "" H 1000 4500 50  0001 C CNN
	1    1000 4500
	1    0    0    -1  
$EndComp
$Comp
L Device:C C3
U 1 1 5C5C0127
P 4250 5100
F 0 "C3" H 3950 5150 50  0000 L CNN
F 1 "1uF" H 3950 5000 50  0000 L CNN
F 2 "Capacitor_SMD:C_0603_1608Metric" H 4288 4950 50  0001 C CNN
F 3 "~" H 4250 5100 50  0001 C CNN
	1    4250 5100
	1    0    0    -1  
$EndComp
Wire Wire Line
	4250 4950 4450 4950
$Comp
L power:GND #PWR0112
U 1 1 5C5C0278
P 4250 5250
F 0 "#PWR0112" H 4250 5000 50  0001 C CNN
F 1 "GND" H 4255 5077 50  0000 C CNN
F 2 "" H 4250 5250 50  0001 C CNN
F 3 "" H 4250 5250 50  0001 C CNN
	1    4250 5250
	1    0    0    -1  
$EndComp
$Comp
L Device:R R3
U 1 1 5C5C088F
P 1750 3950
F 0 "R3" V 1650 3950 50  0000 C CNN
F 1 "22" V 1750 3950 50  0000 C CNN
F 2 "Resistor_SMD:R_0603_1608Metric" V 1680 3950 50  0001 C CNN
F 3 "~" H 1750 3950 50  0001 C CNN
	1    1750 3950
	0    1    1    0   
$EndComp
$Comp
L Device:R R4
U 1 1 5C5C0944
P 1750 4200
F 0 "R4" V 1850 4200 50  0000 C CNN
F 1 "22" V 1750 4200 50  0000 C CNN
F 2 "Resistor_SMD:R_0603_1608Metric" V 1680 4200 50  0001 C CNN
F 3 "~" H 1750 4200 50  0001 C CNN
	1    1750 4200
	0    -1   -1   0   
$EndComp
Wire Wire Line
	1300 4100 1450 4100
Wire Wire Line
	1450 4100 1450 3950
Wire Wire Line
	1450 3950 1600 3950
Wire Wire Line
	1300 4200 1600 4200
Text Label 1900 3950 0    50   ~ 0
D+
Text Label 1900 4200 0    50   ~ 0
D-
Text Label 4450 4650 2    50   ~ 0
D+
Text Label 4450 4750 2    50   ~ 0
D-
Text Label 1300 3900 0    50   ~ 0
VUSB
Text Label 4450 4450 2    50   ~ 0
VUSB
$Comp
L Regulator_Linear:MIC5219YM5 5V_REG1
U 1 1 5C5C2CB8
P 2950 1100
F 0 "5V_REG1" H 2950 1442 50  0000 C CNN
F 1 "MIC5219YM5" H 2950 1351 50  0000 C CNN
F 2 "Package_TO_SOT_SMD:SOT-23-5" H 2950 1425 50  0001 C CNN
F 3 "http://ww1.microchip.com/downloads/en/DeviceDoc/MIC5219-500mA-Peak-Output-LDO-Regulator-DS20006021A.pdf" H 2950 1100 50  0001 C CNN
	1    2950 1100
	1    0    0    -1  
$EndComp
Wire Wire Line
	2350 1000 2650 1000
Wire Wire Line
	3250 1050 3250 1000
$Comp
L Device:R R1
U 1 1 5C5C36A5
P 4400 1250
F 0 "R1" V 4500 1250 50  0000 C CNN
F 1 "1K" V 4400 1250 50  0000 C CNN
F 2 "Resistor_SMD:R_0201_0603Metric" V 4330 1250 50  0001 C CNN
F 3 "~" H 4400 1250 50  0001 C CNN
	1    4400 1250
	-1   0    0    1   
$EndComp
$Comp
L power:GND #PWR0113
U 1 1 5C5C39C9
P 4400 1700
F 0 "#PWR0113" H 4400 1450 50  0001 C CNN
F 1 "GND" H 4405 1527 50  0000 C CNN
F 2 "" H 4400 1700 50  0001 C CNN
F 3 "" H 4400 1700 50  0001 C CNN
	1    4400 1700
	1    0    0    -1  
$EndComp
Wire Wire Line
	4400 1100 3750 1100
Wire Wire Line
	3750 1100 3750 1050
Connection ~ 3750 1050
Wire Wire Line
	3750 1050 3950 1050
$Comp
L Device:LED D1
U 1 1 5C5C3E88
P 4400 1550
F 0 "D1" V 4450 1750 50  0000 R CNN
F 1 "LED" V 4350 1750 50  0000 R CNN
F 2 "LED_SMD:LED_0603_1608Metric" H 4400 1550 50  0001 C CNN
F 3 "~" H 4400 1550 50  0001 C CNN
	1    4400 1550
	0    -1   -1   0   
$EndComp
Wire Notes Line
	1900 2000 4600 2000
Wire Notes Line
	1900 600  4600 600 
$Comp
L power:GND #PWR0114
U 1 1 5C5C4B32
P 3250 1650
F 0 "#PWR0114" H 3250 1400 50  0001 C CNN
F 1 "GND" H 3255 1477 50  0000 C CNN
F 2 "" H 3250 1650 50  0001 C CNN
F 3 "" H 3250 1650 50  0001 C CNN
	1    3250 1650
	1    0    0    -1  
$EndComp
$Comp
L Device:C C1
U 1 1 5C5C4B38
P 3250 1500
F 0 "C1" H 3365 1546 50  0000 L CNN
F 1 "470pF" H 3365 1455 50  0000 L CNN
F 2 "Capacitor_SMD:C_0603_1608Metric" H 3288 1350 50  0001 C CNN
F 3 "~" H 3250 1500 50  0001 C CNN
	1    3250 1500
	1    0    0    -1  
$EndComp
Wire Wire Line
	3250 1350 3250 1200
Connection ~ 3750 1100
Wire Wire Line
	3750 1100 3750 1350
Wire Wire Line
	3250 1050 3750 1050
Wire Wire Line
	2350 1000 2350 1100
Wire Wire Line
	2350 1100 2650 1100
Connection ~ 2350 1100
Wire Wire Line
	2350 1100 2350 1350
$Comp
L Device:R R0
U 1 1 5C5C6D45
P 4300 3500
F 0 "R0" V 4400 3500 50  0000 C CNN
F 1 "1K" V 4300 3500 50  0000 C CNN
F 2 "Resistor_SMD:R_0201_0603Metric" V 4230 3500 50  0001 C CNN
F 3 "~" H 4300 3500 50  0001 C CNN
	1    4300 3500
	-1   0    0    1   
$EndComp
Wire Wire Line
	4300 3650 4450 3650
Text Label 3950 1050 0    50   ~ 0
VCC
Text Label 5050 3350 1    50   ~ 0
VCC
Text Label 4300 3350 1    50   ~ 0
VCC
Wire Notes Line
	5000 600  8400 600 
Wire Notes Line
	8400 600  8400 2500
Wire Notes Line
	8400 2500 5000 2500
Wire Notes Line
	5000 2500 5000 600 
Text Notes 8250 750  2    50   ~ 0
Solenoid\n
Wire Notes Line
	11000 600  8800 600 
Text Notes 10900 750  2    50   ~ 0
Pressure Sensors
$Comp
L Device:Ferrite_Bead FB1
U 1 1 5C5CE769
P 2350 5400
F 0 "FB1" H 2487 5446 50  0000 L CNN
F 1 "Ferrite_Bead" H 2487 5355 50  0000 L CNN
F 2 "Resistor_SMD:R_0402_1005Metric" V 2280 5400 50  0001 C CNN
F 3 "~" H 2350 5400 50  0001 C CNN
	1    2350 5400
	1    0    0    -1  
$EndComp
$Comp
L Device:C C5
U 1 1 5C5CE8A2
P 2350 5700
F 0 "C5" H 2465 5746 50  0000 L CNN
F 1 "100nF" H 2465 5655 50  0000 L CNN
F 2 "Capacitor_SMD:C_0603_1608Metric" H 2388 5550 50  0001 C CNN
F 3 "~" H 2350 5700 50  0001 C CNN
	1    2350 5700
	1    0    0    -1  
$EndComp
Text Label 5150 3350 1    50   ~ 0
AVCC
Text Label 2000 5550 2    50   ~ 0
AVCC
Wire Wire Line
	2000 5550 2350 5550
Connection ~ 2350 5550
$Comp
L power:GND #PWR0115
U 1 1 5C5CF1F1
P 2350 5850
F 0 "#PWR0115" H 2350 5600 50  0001 C CNN
F 1 "GND" H 2355 5677 50  0000 C CNN
F 2 "" H 2350 5850 50  0001 C CNN
F 3 "" H 2350 5850 50  0001 C CNN
	1    2350 5850
	1    0    0    -1  
$EndComp
Text Label 2350 5250 1    50   ~ 0
VCC
Text Label 4950 3350 1    50   ~ 0
VCC
Text Label 850  5500 1    50   ~ 0
VUSB
$Comp
L Device:C C4
U 1 1 5C5D01EA
P 850 5650
F 0 "C4" H 965 5696 50  0000 L CNN
F 1 "10uF" H 965 5605 50  0000 L CNN
F 2 "Capacitor_SMD:C_0603_1608Metric" H 888 5500 50  0001 C CNN
F 3 "~" H 850 5650 50  0001 C CNN
	1    850  5650
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR0116
U 1 1 5C5D0268
P 850 5800
F 0 "#PWR0116" H 850 5550 50  0001 C CNN
F 1 "GND" H 855 5627 50  0000 C CNN
F 2 "" H 850 5800 50  0001 C CNN
F 3 "" H 850 5800 50  0001 C CNN
	1    850  5800
	1    0    0    -1  
$EndComp
Text Label 4450 4250 2    50   ~ 0
AREF
Text Label 1400 5500 1    50   ~ 0
AREF
$Comp
L Device:C C6
U 1 1 5C5D0EFC
P 1400 5650
F 0 "C6" H 1515 5696 50  0000 L CNN
F 1 "1uF" H 1515 5605 50  0000 L CNN
F 2 "Capacitor_SMD:C_0603_1608Metric" H 1438 5500 50  0001 C CNN
F 3 "~" H 1400 5650 50  0001 C CNN
	1    1400 5650
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR0117
U 1 1 5C5D0F8B
P 1400 5800
F 0 "#PWR0117" H 1400 5550 50  0001 C CNN
F 1 "GND" H 1405 5627 50  0000 C CNN
F 2 "" H 1400 5800 50  0001 C CNN
F 3 "" H 1400 5800 50  0001 C CNN
	1    1400 5800
	1    0    0    -1  
$EndComp
$Comp
L Device:Crystal Y1
U 1 1 5C5D218C
P 3900 3950
F 0 "Y1" V 3854 4081 50  0000 L CNN
F 1 "16MHz" V 3945 4081 50  0000 L CNN
F 2 "SpiroParts:Crystal_SMD_SeikoEpson_TSX3225-4Pin_3.2x2.5mm" H 3900 3950 50  0001 C CNN
F 3 "~" H 3900 3950 50  0001 C CNN
	1    3900 3950
	0    1    1    0   
$EndComp
Wire Wire Line
	4350 3800 4350 3850
Wire Wire Line
	4350 3850 4450 3850
Wire Wire Line
	4350 4100 4350 4050
Wire Wire Line
	4350 4050 4450 4050
$Comp
L Device:C C7
U 1 1 5C5D444C
P 3500 3800
F 0 "C7" V 3752 3800 50  0000 C CNN
F 1 "22pF" V 3661 3800 50  0000 C CNN
F 2 "Capacitor_SMD:C_0603_1608Metric" H 3538 3650 50  0001 C CNN
F 3 "~" H 3500 3800 50  0001 C CNN
	1    3500 3800
	0    -1   -1   0   
$EndComp
$Comp
L Device:C C8
U 1 1 5C5D454B
P 3500 4100
F 0 "C8" V 3650 4100 50  0000 C CNN
F 1 "22pF" V 3750 4100 50  0000 C CNN
F 2 "Capacitor_SMD:C_0603_1608Metric" H 3538 3950 50  0001 C CNN
F 3 "~" H 3500 4100 50  0001 C CNN
	1    3500 4100
	0    1    1    0   
$EndComp
Wire Wire Line
	3350 4100 3350 3950
Wire Wire Line
	3650 3800 3900 3800
Wire Wire Line
	3900 3800 4350 3800
Connection ~ 3900 3800
Wire Wire Line
	4350 4100 3900 4100
Connection ~ 3900 4100
Wire Wire Line
	3900 4100 3650 4100
$Comp
L power:GND #PWR0118
U 1 1 5C5D5617
P 3350 3950
F 0 "#PWR0118" H 3350 3700 50  0001 C CNN
F 1 "GND" V 3355 3822 50  0000 R CNN
F 2 "" H 3350 3950 50  0001 C CNN
F 3 "" H 3350 3950 50  0001 C CNN
	1    3350 3950
	0    1    1    0   
$EndComp
Connection ~ 3350 3950
Wire Wire Line
	3350 3950 3350 3800
Wire Notes Line
	600  3000 6550 3000
Wire Notes Line
	6550 3000 6550 7650
Wire Notes Line
	6550 7650 600  7650
Text Notes 6400 3150 2    50   ~ 0
Microcontroller
Wire Wire Line
	900  4500 1000 4500
Connection ~ 1000 4500
$Comp
L Connector:TestPoint TP1
U 1 1 5C5EE45F
P 6900 3850
F 0 "TP1" H 6958 3970 50  0000 L CNN
F 1 "TestPoint" H 6958 3879 50  0000 L CNN
F 2 "SpiroParts:testPoint-2.2x1.4" H 7100 3850 50  0001 C CNN
F 3 "~" H 7100 3850 50  0001 C CNN
	1    6900 3850
	1    0    0    -1  
$EndComp
$Comp
L Connector:TestPoint TP2
U 1 1 5C5EE560
P 7450 3850
F 0 "TP2" H 7508 3970 50  0000 L CNN
F 1 "TestPoint" H 7508 3879 50  0000 L CNN
F 2 "SpiroParts:testPoint-2.2x1.4" H 7650 3850 50  0001 C CNN
F 3 "~" H 7650 3850 50  0001 C CNN
	1    7450 3850
	1    0    0    -1  
$EndComp
$Comp
L Connector:TestPoint TP3
U 1 1 5C5EE608
P 8000 3850
F 0 "TP3" H 8058 3970 50  0000 L CNN
F 1 "TestPoint" H 8058 3879 50  0000 L CNN
F 2 "SpiroParts:testPoint-2.2x1.4" H 8200 3850 50  0001 C CNN
F 3 "~" H 8200 3850 50  0001 C CNN
	1    8000 3850
	1    0    0    -1  
$EndComp
Text Label 8000 3850 3    50   ~ 0
PressureSensorLow2
Text Label 7450 3850 3    50   ~ 0
PressureSensorLow1
Text Label 6900 3850 3    50   ~ 0
PressureSensorHigh
$Comp
L Connector:TestPoint TP4
U 1 1 5C5F39CA
P 6900 4950
F 0 "TP4" H 6958 5070 50  0000 L CNN
F 1 "TestPoint" H 6958 4979 50  0000 L CNN
F 2 "SpiroParts:testPoint-2.2x1.4" H 7100 4950 50  0001 C CNN
F 3 "~" H 7100 4950 50  0001 C CNN
	1    6900 4950
	1    0    0    -1  
$EndComp
$Comp
L Connector:TestPoint TP5
U 1 1 5C5F3A76
P 7500 4950
F 0 "TP5" H 7558 5070 50  0000 L CNN
F 1 "TestPoint" H 7558 4979 50  0000 L CNN
F 2 "SpiroParts:testPoint-2.2x1.4" H 7700 4950 50  0001 C CNN
F 3 "~" H 7700 4950 50  0001 C CNN
	1    7500 4950
	1    0    0    -1  
$EndComp
$Comp
L Connector:TestPoint TP6
U 1 1 5C5F3B1D
P 8050 4950
F 0 "TP6" H 8108 5070 50  0000 L CNN
F 1 "TestPoint" H 8108 4979 50  0000 L CNN
F 2 "SpiroParts:testPoint-2.2x1.4" H 8250 4950 50  0001 C CNN
F 3 "~" H 8250 4950 50  0001 C CNN
	1    8050 4950
	1    0    0    -1  
$EndComp
Text Label 6900 4950 3    50   ~ 0
VUSB
$Comp
L Connector:TestPoint TP0
U 1 1 5C5F47A8
P 6900 3250
F 0 "TP0" H 6958 3370 50  0000 L CNN
F 1 "TestPoint" H 6958 3279 50  0000 L CNN
F 2 "SpiroParts:testPoint-2.2x1.4" H 7100 3250 50  0001 C CNN
F 3 "~" H 7100 3250 50  0001 C CNN
	1    6900 3250
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR0119
U 1 1 5C5F48DC
P 6900 3250
F 0 "#PWR0119" H 6900 3000 50  0001 C CNN
F 1 "GND" H 6905 3077 50  0000 C CNN
F 2 "" H 6900 3250 50  0001 C CNN
F 3 "" H 6900 3250 50  0001 C CNN
	1    6900 3250
	1    0    0    -1  
$EndComp
Text Label 7500 4950 3    50   ~ 0
VCC
Text Label 8050 4950 3    50   ~ 0
solenoidToggle
Wire Notes Line
	8800 600  8800 4700
Wire Notes Line
	8800 4700 11000 4700
Wire Notes Line
	11000 4700 11000 600 
$Comp
L Connector:TestPoint TP7
U 1 1 5C5FAC4B
P 6900 5800
F 0 "TP7" H 6958 5920 50  0000 L CNN
F 1 "TestPoint" H 6958 5829 50  0000 L CNN
F 2 "SpiroParts:testPoint-2.2x1.4" H 7100 5800 50  0001 C CNN
F 3 "~" H 7100 5800 50  0001 C CNN
	1    6900 5800
	1    0    0    -1  
$EndComp
Text Label 6900 5800 3    50   ~ 0
Vin
$Comp
L Mechanical:MountingHole H2
U 1 1 5C5FB11D
P 10300 5100
F 0 "H2" H 10400 5146 50  0000 L CNN
F 1 "MountingHole" H 10400 5055 50  0000 L CNN
F 2 "MountingHole:MountingHole_3.2mm_M3" H 10300 5100 50  0001 C CNN
F 3 "~" H 10300 5100 50  0001 C CNN
	1    10300 5100
	1    0    0    -1  
$EndComp
$Comp
L Mechanical:MountingHole H1
U 1 1 5C5FB344
P 9500 5100
F 0 "H1" H 9600 5146 50  0000 L CNN
F 1 "MountingHole" H 9600 5055 50  0000 L CNN
F 2 "MountingHole:MountingHole_3.2mm_M3" H 9500 5100 50  0001 C CNN
F 3 "~" H 9500 5100 50  0001 C CNN
	1    9500 5100
	1    0    0    -1  
$EndComp
$Comp
L Mechanical:MountingHole H3
U 1 1 5C5FB3F8
P 9500 5400
F 0 "H3" H 9600 5446 50  0000 L CNN
F 1 "MountingHole" H 9600 5355 50  0000 L CNN
F 2 "MountingHole:MountingHole_3.2mm_M3" H 9500 5400 50  0001 C CNN
F 3 "~" H 9500 5400 50  0001 C CNN
	1    9500 5400
	1    0    0    -1  
$EndComp
$Comp
L Mechanical:MountingHole H4
U 1 1 5C5FB4A9
P 10300 5400
F 0 "H4" H 10400 5446 50  0000 L CNN
F 1 "MountingHole" H 10400 5355 50  0000 L CNN
F 2 "MountingHole:MountingHole_3.2mm_M3" H 10300 5400 50  0001 C CNN
F 3 "~" H 10300 5400 50  0001 C CNN
	1    10300 5400
	1    0    0    -1  
$EndComp
Wire Wire Line
	1300 1550 1300 1650
Connection ~ 1300 1650
$Comp
L Connector:TestPoint TP9
U 1 1 5C5FF6A3
P 8050 5800
F 0 "TP9" H 8108 5920 50  0000 L CNN
F 1 "TestPoint" H 8108 5829 50  0000 L CNN
F 2 "SpiroParts:testPoint-2.2x1.4" H 8250 5800 50  0001 C CNN
F 3 "~" H 8250 5800 50  0001 C CNN
	1    8050 5800
	1    0    0    -1  
$EndComp
Text Label 8050 5800 3    50   ~ 0
solenoidGND
Text Label 900  6500 1    50   ~ 0
VCC
$Comp
L Device:C C9
U 1 1 5C63F564
P 900 6650
F 0 "C9" H 1015 6696 50  0000 L CNN
F 1 "0.1uF" H 1015 6605 50  0000 L CNN
F 2 "Capacitor_SMD:C_0603_1608Metric" H 938 6500 50  0001 C CNN
F 3 "~" H 900 6650 50  0001 C CNN
	1    900  6650
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR0120
U 1 1 5C63F56B
P 900 6800
F 0 "#PWR0120" H 900 6550 50  0001 C CNN
F 1 "GND" H 905 6627 50  0000 C CNN
F 2 "" H 900 6800 50  0001 C CNN
F 3 "" H 900 6800 50  0001 C CNN
	1    900  6800
	1    0    0    -1  
$EndComp
Text Label 1850 6450 1    50   ~ 0
VCC
$Comp
L Device:C C10
U 1 1 5C647D5B
P 1850 6600
F 0 "C10" H 1965 6646 50  0000 L CNN
F 1 "0.1uF" H 1965 6555 50  0000 L CNN
F 2 "Capacitor_SMD:C_0603_1608Metric" H 1888 6450 50  0001 C CNN
F 3 "~" H 1850 6600 50  0001 C CNN
	1    1850 6600
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR0121
U 1 1 5C647D62
P 1850 6750
F 0 "#PWR0121" H 1850 6500 50  0001 C CNN
F 1 "GND" H 1855 6577 50  0000 C CNN
F 2 "" H 1850 6750 50  0001 C CNN
F 3 "" H 1850 6750 50  0001 C CNN
	1    1850 6750
	1    0    0    -1  
$EndComp
$Comp
L Transistor_FET:2N7002 NMOS1
U 1 1 5C64C233
P 7400 1350
F 0 "NMOS1" H 7605 1396 50  0000 L CNN
F 1 "2N7002" H 7605 1305 50  0000 L CNN
F 2 "Package_TO_SOT_SMD:SOT-23" H 7600 1275 50  0001 L CIN
F 3 "https://www.fairchildsemi.com/datasheets/2N/2N7002.pdf" H 7400 1350 50  0001 L CNN
	1    7400 1350
	1    0    0    -1  
$EndComp
Text Label 2550 6500 1    50   ~ 0
VCC
$Comp
L Device:C C11
U 1 1 5C655CD2
P 2550 6650
F 0 "C11" H 2665 6696 50  0000 L CNN
F 1 "0.1uF" H 2665 6605 50  0000 L CNN
F 2 "Capacitor_SMD:C_0603_1608Metric" H 2588 6500 50  0001 C CNN
F 3 "~" H 2550 6650 50  0001 C CNN
	1    2550 6650
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR0122
U 1 1 5C655CD9
P 2550 6800
F 0 "#PWR0122" H 2550 6550 50  0001 C CNN
F 1 "GND" H 2555 6627 50  0000 C CNN
F 2 "" H 2550 6800 50  0001 C CNN
F 3 "" H 2550 6800 50  0001 C CNN
	1    2550 6800
	1    0    0    -1  
$EndComp
Text Label 3200 6500 1    50   ~ 0
AVCC
$Comp
L Device:C C12
U 1 1 5C6563FF
P 3200 6650
F 0 "C12" H 3315 6696 50  0000 L CNN
F 1 "0.1uF" H 3315 6605 50  0000 L CNN
F 2 "Capacitor_SMD:C_0603_1608Metric" H 3238 6500 50  0001 C CNN
F 3 "~" H 3200 6650 50  0001 C CNN
	1    3200 6650
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR0123
U 1 1 5C656406
P 3200 6800
F 0 "#PWR0123" H 3200 6550 50  0001 C CNN
F 1 "GND" H 3205 6627 50  0000 C CNN
F 2 "" H 3200 6800 50  0001 C CNN
F 3 "" H 3200 6800 50  0001 C CNN
	1    3200 6800
	1    0    0    -1  
$EndComp
Text Label 1400 6900 1    50   ~ 0
VCC
$Comp
L Device:C C13
U 1 1 5C65725B
P 1400 7050
F 0 "C13" H 1515 7096 50  0000 L CNN
F 1 "0.1uF" H 1515 7005 50  0000 L CNN
F 2 "Capacitor_SMD:C_0603_1608Metric" H 1438 6900 50  0001 C CNN
F 3 "~" H 1400 7050 50  0001 C CNN
	1    1400 7050
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR0124
U 1 1 5C657262
P 1400 7200
F 0 "#PWR0124" H 1400 6950 50  0001 C CNN
F 1 "GND" H 1405 7027 50  0000 C CNN
F 2 "" H 1400 7200 50  0001 C CNN
F 3 "" H 1400 7200 50  0001 C CNN
	1    1400 7200
	1    0    0    -1  
$EndComp
$EndSCHEMATC