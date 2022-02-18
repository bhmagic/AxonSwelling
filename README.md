# AxonSwelling
hybrid 1D-2D model for AP passing axon swelling 

The code provided here is the exact code we used in this article. The code is programed and tested in Matlab 2014b. The main code is named “Hybrid.m”, and there are four other required functions files (initBCimplicit_1.m, initBCimplicit_head.m, initBCimplicit_tail.m, loopBC_1.m) in the same zip file. All the parameters and variable initialization can be setup at the beginning of the code. 

Running the code will automatically plot the model layout, the meshing, and a real-time monitoring of all the status, as showing in Figure S2. At the end of simulation, the history of the GUI will be documented as a GIF file for ease of access and all the physical variables and other intermediate variables are saved in two files (data.mat and all.mat) for later uses. 


Wu, Y. T., Gilpin, K., & Adnan, A. (2020). Effects of focal axonal swelling level on the action potential signal transmission. Journal of computational neuroscience, 48(3), 253-263.
