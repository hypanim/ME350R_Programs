#include <Arduino.h>

int relay1_pin = 5;
int relay2_pin = 6;
int delay_time_ms = 25000;
int cooldown_ms = 10000;

void setup() {
  Serial.begin(115200);
  pinMode(relay1_pin, OUTPUT);
  pinMode(relay2_pin, OUTPUT);
}

void loop() {
  Serial.println("running");
  digitalWrite(relay1_pin, HIGH);
  digitalWrite(relay2_pin, LOW);
  delay(delay_time_ms);
  digitalWrite(relay1_pin, LOW);
  digitalWrite(relay2_pin, LOW);
  delay(cooldown_ms);
  digitalWrite(relay1_pin, LOW);
  digitalWrite(relay2_pin, HIGH);
  delay(delay_time_ms);
  digitalWrite(relay1_pin, LOW);
  digitalWrite(relay2_pin, LOW);
  delay(cooldown_ms);
}
