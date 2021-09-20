# Quantum-dot-laser-applications
This project is to apply quantum dot laser in microwave signal processing. 

1. Microwave photonics filter simulation - include two classes: 
        
         a. comb generatio
        
         b. MWP filter response
2. Comb shaping using waveshaper - function includes:

         a. Flatten comb

         b. Gaussian apodization

         c. Sampling

         d. Sinc for flat-top filter response

         f. Sinc^2 for triangular filter response

         e. Sinc*cos for flat-top bandpass filter

         g. Sinc^2*cos for triangular bandpass filter
   
 3. Read ESA filter respones data 

         Sweep range: 0-2.5GHz/2.5-5GHz/5-7.5GHz; combine these three data in code
 6. Compare sim and exp
   
         Step1: Read ESA data 

         Step2: Run the code - include: 

               a. Load measured OFC data and ESA data

               b. Simulate using measured OFC and compare with the ESA data
   
