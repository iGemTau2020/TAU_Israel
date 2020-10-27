Contained within this directory are scripts used for the optimizer, and the methylation database used within the software.


The code is still presented, and many internal functions from it can be scavenged for future projects. For any help, feel free to contact us at igem.tau.2019@gmail.com.


Staubility EFM Optimizer Code
------------------------------
This folder contains the source code of our product, the Staubility EFM Optimizer.

Our code combines many frameworks, some of whom only work on certain operating systems, and some of whom require highly nontrivial installations and dependencies. Thus, running our code as a script is extremely tricky. This is part of the reason we built the GUI framework, to avoid having our users go through this process, and compiled it into a software.exe, so our users can run it without creating and installing the different dependencies.

Our software.exe was compiled from the code in this folder which is specified in the next section.

Folder Contents
----------------
 - EFM_New.py - This code contains the different models of the project (SSR Detector, RMD detector Methylation Detector, sequence_optimizer) inside one file with few alignments for running inside a complex program, and a main function which activates the models according to the user inputs and the program step.
 
 - Gui.py -  The Gui code is our software interface code. The user interface was created with the unique framework we built for this project, based on the pyqt library. The Gui_code file contains the customer side functions which test and modify the different inputs, the user interfaces structured code and our framework.

