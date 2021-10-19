# DSP_SG

1. Generate OFDM 20MHz,
![image](https://user-images.githubusercontent.com/87049112/135781213-354d0766-5c56-4638-8c2a-500a2aec7845.png)

2. Apply Channel filter and gain, evm: 0.12, delay: 57
![image](https://user-images.githubusercontent.com/87049112/136729273-34726f42-5641-4ee9-83ee-7f932017ecf5.png)

3. Generate HB and SRC FIR,

database FIR:

![image](https://user-images.githubusercontent.com/87049112/136647623-587d3d74-84a7-4c42-83c7-2e6901f395e6.png)

generate FIR(red) comparsion:
![image](https://user-images.githubusercontent.com/87049112/136727516-18e255ed-845a-41a1-b501-2bd88b8ee11f.png)

4. Interpolation to fsOut=122.88e6
![image](https://user-images.githubusercontent.com/87049112/136744489-a16c2439-b897-4211-86c2-4399730ffac7.png)

5. Export signal to Fixed point
![image](https://user-images.githubusercontent.com/87049112/136744601-b7afbc5c-ada7-497f-915d-d2b48213725f.png)

6. Synchronize and check EVM=1.64%

7. Add NCO = [-20e6 20e6] and combine carriers
![image](https://user-images.githubusercontent.com/87049112/137058117-59fa3db5-681c-4545-bf90-5146a54be0a0.png)

8. PAR comparsion
![image](https://user-images.githubusercontent.com/87049112/137826661-bb43dea7-641b-4d9a-b5d5-5d6dbe438849.png)

