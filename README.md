
# Digit Recognition from Speech

## Overview
This project focuses on recognizing spoken digits by using Linear Predictive Coding (LPC) coefficients and Hidden Markov Models. The system processes audio files, classifies the spoken digits, and stores results in designated output folders.

## Getting Started

### Prerequisites
- **Microsoft Visual Studio** for building and running the project.
- **C/C++** programming knowledge.

### Steps to Open and Run the Project

1. **Open the Solution File:**
   - Navigate to the "234101006_digit_recognition" folder.
   - Open the `234101006_digit_recognition.sln` file in Microsoft Visual Studio.

2. **Build the Project:**
   - Press `F5` to build the project. This will compile the `234101006_digit_recognition.cpp` file.

3. **Compile the Code:**
   - After building, press `F7` to compile the code.

4. **Testing and Results:**
   - The output window will display testing results based on 10 files of utterances for each digit.
   - Lambda models for each digit sequence are saved in text files located in the `output` folder, under respective digit number folders.
   - Observation sequences and average models for each digit are also stored in the `output` folder.
   - Testing results are saved in text files located in the `results` folder.

## Files and Folders

- **234101006_digit_recognition.sln**: Solution file for Microsoft Visual Studio.
- **234101006_digit_recognition.cpp**: Source code for the digit recognition system.
- **output/**: Contains lambda models, observation sequences, and average models.
- **results/**: Contains the testing results.

