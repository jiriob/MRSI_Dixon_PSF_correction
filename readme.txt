You will need to compile the function ba_interp3.cpp available here: http://www.mathworks.at/matlabcentral/fileexchange/21702-3d-volume-interpolation-with-ba-interp3--fast-interp3-replacement

In my experience, it's not easy to have a functioning c++ compiler working in Matlab. You have to have everything uptodate (the newest Matlab and OS), or you can try to find and install compiler compatible with older versions of Matlab, if you manage to and don't have the option to update Matlab. 

The installation process is done in Matlab using this command:

mex '/path_to_the_ba_interp3_folder/ba_interp3.cpp'

Be sure you don't have the "ba_interp3.m" file in your Matlab path,) because than, Matlab will try to use the .m file instead of the compiled one (which is saved in the defauld Matlab folder in Documents).


After the compilation you have to divide your data into three folders:
- Dixon images are in two of them: water dicoms are in a folder called "Dixon_W" and lipid dicoms in "Dixon_F"
- the 3D MRSI dicom from Siemens is saved in a folder named "Spec"

You have an option to calculate the SNR in your MRSI using the function SNR_lenk.m .



