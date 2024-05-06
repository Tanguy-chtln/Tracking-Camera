
// #include <Arduino.h>
// const int WRITE_PIN = 12;
// const int READ_PIN = 13;

// unsigned long time;
// unsigned long lastTime;
// bool doWrite = false;
// String text("");
// const int MAX_LENGTH = 20;
// char incoming_Byte = '\r';
// bool is_end_of_line(char c) { return c == '\r' or c == '\n'; }
// bool read_line()
// {
//   if (is_end_of_line(incoming_Byte) or text.length() == MAX_LENGTH)
//   {
//     text = "";
//     incoming_Byte = 'a';
//   }
//   while (Serial.available() and text.length() < MAX_LENGTH)
//   {
//     incoming_Byte = Serial.read();
//     if (is_end_of_line(incoming_Byte))
//     {
//       break;
//     }
//     text += incoming_Byte;
//   }
//   return is_end_of_line(incoming_Byte) or text.length() == MAX_LENGTH;
// }

// void setup()
// {
//     time = 0;
//     lastTime = 0;
//     Serial.begin(19200);
//     pinMode(WRITE_PIN, OUTPUT);
//     pinMode(READ_PIN, OUTPUT);
//     digitalWrite(WRITE_PIN, LOW);
//     digitalWrite(READ_PIN, LOW);
// }


// void loop()
// {
//     time = millis();
//     if (time - lastTime > 1000) {
//         if (doWrite) {
//             digitalWrite(READ_PIN, LOW);
//             digitalWrite(WRITE_PIN, HIGH);
//             Serial.println(text);
//             lastTime = time;
//             doWrite = !doWrite;
//         } else {
//             digitalWrite(WRITE_PIN, LOW);
//             digitalWrite(READ_PIN, HIGH);
//             if (read_line()) {
//                 lastTime = time;
//                 doWrite = !doWrite;
//             }

//         }
//     }

// }  