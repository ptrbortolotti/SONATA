Introduction
------------------------------------------------------------------------------------------------------
This installation package is for VABS, an efficient high-fidelity cross-sectional analysis tool
for design and analysis of composite beams. More details can be found in VABSManual.pdf. 

VABS is a commercial code. You need to obtain a license (trial, academic, commercial) from 
AnalySwift, LLC (http://analyswift.com/). 
------------------------------------------------------------------------------------------------------


Get Started
------------------------------------------------------------------------------------------------------
Double clicking the installation executable and follow the setup wizard will install SwiftComp. 
After installation, you need to do the following three steps to run the code.  
1. Add the installation folder to the path of the system or user. This can be done through editing
    the Environmental Variable: PATH. This enables the system to locate VABSIII so that you 
    can execute the code in any folder. 
2. Type VABSIII, it will output: Something wrong with the license file:  Can't connect to license 
    server "xxxx", the name inside the quotes is your computer (host) name. Use this to request a 
    trial license from AnalySwift. Put the license file inside the installation folder.
3. You are now to ready to use VABS. To verify whether you have successfully installed VABS, 
    use command line window to enter the example folder and type "VABSIII isorect.dat". 
    If everything is configured right, VABS will report success of the run.  
------------------------------------------------------------------------------------------------------


Prepare Input Files
------------------------------------------------------------------------------------------------------
VABS is a command line code, and you can prepare your own input file according 
to VABSManual.pdf. Input files for several simple cross sections along with a ppt 
file showing the cross sections can be found in folder "examples". You can also 
use PreVABS, an open-source, parametric pre- and post-processor for VABS. 
PreVABS can be downloaded from https://cdmhub.org/resources/1597. You 
can also check with AnalySwift for interfaces with other software programs.  
------------------------------------------------------------------------------------------------------


Official Licenses
------------------------------------------------------------------------------------------------------
After the trial period, you can request either a node-locked license or a floating license. 
1. To request a node-locked license, all you need to do is to provide your username and put the 
    received license file inside the installation folder.
2. To request a floating license, you need to provide the computer (host) name and the mac address. Put 
    the received license file inside the installation folder and execute sglmserver license_file_name to start
    the license manager first. It is recommended to keep the license manager running as you need to
    validate the license every time you run VABSIII.
3. If you want to manage floating licenses using a computer (a server) different from the computers 
   (clients) on which you run the code, you need to do the following: 
   a. Install the code both on the server and the client.
   b. Provide the server's name and mac address to request a valid license file from AnalySwift. 
   c. Put the license file inside the installation folder on the server.  
   d. Open port 29750 on the server to allow client machines to talk to the server. 
   e. execute sglmserver license_file_name to start the license server. 
   f.  Add an environmental variable "SG_LICENSE_FILE" with value "29750@servername" on client
       machines. 
   g. Verify the configuration of license manager by running VABSIII on client machines. 
------------------------------------------------------------------------------------------------------


Uninstallation
------------------------------------------------------------------------------------------------------
If you decide to remove VABS, please double click uninstall.exe in the installation folder. After 
uninstallation, you will be directed to the website of AnalySwift.com to provide some feedbacks. 
------------------------------------------------------------------------------------------------------


EULA Acceptance
------------------------------------------------------------------------------------------------------
You are directed to the End User License Agreement (the "EULA"), which is included in the download 
as "EULA.txt". By downloading, executing, or running this software, I represent that (i) I have 
read AnalySwift's End User License Agreement (the "EULA") in its entirety; (ii) I, or the entity 
on whose behalf I am acting, agree to be bound by the terms and conditions of the EULA; and (iii) 
if I am acting on behalf of an entity, I am authorized to accept the EULA on behalf of such entity.
------------------------------------------------------------------------------------------------------


Tech Support
------------------------------------------------------------------------------------------------------
Please signup our newsletter at http://analyswift.com/ to receive updates and news regarding VABS. 
Please join Dr. Wenbin Yu's Group in the Cloud at https://cdmhub.org/groups/yugroup and post your
not-proprietary questions on the forum within the group (https://cdmhub.org/groups/yugroup/forum). 
------------------------------------------------------------------------------------------------------