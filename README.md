Welcome to QuantumSimLite application made by Maxime Vinteler (@maximev1314 on Twitter) ! This document provides essential information to help you navigate and make the most of the functionalities offered by my software.

## Key Features:

1. **Eigenstates Visualization:**
   - Explore eigenstates for any potential in 1D and 2D for a single particle.
   - View detailed visualizations in the application.

2. **Quantum Dynamics:**
   - Conduct quantum dynamics simulations within the application.
   - Experience real-time changes in the system.

3. **Downloadable Resources:**
   - Download eigenstates and dynamic simulations (FFMPEG require for videos).
   - Note: Saving images or videos may cause a temporary freeze in the application.

## Important Information:

- **Compatibility:**
   - This version of the application is designed for Windows, MACOS and Linux operating systems.
   - For Windows users, we also have an .exe version, which boasts the advantage of being easier to install. : https://github.com/MaximeV1314/QuantumSimLite-Windows
   - Access crucial details by hovering over the "blue information icons" in the application.
   - Read and follow the provided information for correct usage of the software.

- **Information Tools:**
   - Access crucial details by hovering over the "blue information icons" in the application.
   - Read and follow the provided information for correct usage of the software.
   - Access to the eigenstates and eigenenergies data in VP folder, .bz2 (eigenstates) and .csv (eigenenergies) files. Note that for 2D, you need to reshape 1D arrays of shape N^2 into an array of 2D shape (NxN).
   - Access to the dynamics data in VP folder --> dynamic folder. Files with "info" as a prefix contain mean value and standard deviation data for position and momentum representations. Note that for 2D, you need to reshape 1D arrays of shape N^2 into an array of 2D shape (NxN).
     
- **Media Folder:**
   - Locate additional resources, including videos and images, in the 'media' folder.
   - Utilize these materials for a richer experience and deeper insights.

 ## System Requirements:

Ensure your system meets the following requirements for optimal performance:
- Operating System: Windows, MACOS, Linux.
- Disk Space: 40 MB of available space

## Installation:

1. Download the repository (Code --> Download ZIP) at https://github.com/MaximeV1314/QuantumSimLite-Win-MacOS-Linux
2. Unzip the file in your directory.
3. Ensure that the following Python libraries are installed on your computer (it is recommended to use the latest version of Python along with these libraries):
    - Matplotlib at https://matplotlib.org/stable/users/installing/index.html
    - Numpy at https://numpy.org/install/
    - Pandas at https://pandas.pydata.org/docs/getting_started/install.html
    - PIL at https://pillow.readthedocs.io/en/stable/installation.html
    - Idlelib wich is already included with python 3.7+ on Windows and MACOS. On Linux, you can install it using the command: sudo apt-get install idle3
4. (Optionnal) FFMPEG is require to download dynamic videos. I recommend following these straightforward tutorials to install it :
   - For Windows : https://phoenixnap.com/kb/ffmpeg-windows
   - For MACOS : https://phoenixnap.com/kb/ffmpeg-mac
   - For Linux : https://phoenixnap.com/kb/install-ffmpeg-ubuntu
5. Run the Python script QuantumSimLite.py (note that the first time you open it, the application may take some time).

## Usage Guidelines:

1. Open the application and explore the various features available.
2. Refer to the information tools for guidance on using specific functionalities.
3. Take advantage of the media folder to access additional resources.

If you encounter any issues or have questions, feel free to reach out at maxime.vinteler@yahoo.fr or @maximev1314 on Twitter.

Special thanks to Paul Rieunier (@LaMatriceCarree and @rotor_mr on Twitter) for creating the application's logo.

Happy Quantum Simulating!

## Version
- 2.0 Add 2D simulatons
- 1.1 Can now put the app in full screen
- 1.0 Release of the application
