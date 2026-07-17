# INSTALL THE APP : 

## 1) Install Docker Desktop
Download docker desktop then open it :

For Windows : 
https://docs.docker.com/desktop/setup/install/windows-install/

For Mac (for older macs choose ): 
https://docs.docker.com/desktop/setup/install/mac-install/


## 2) Pull the docker image from Dockerhub :
Run in the terminal :
> docker pull qtea1/bulktools:latest


# RUN THE APP : 
## FOR WINDOWS (two options) : 
### Click style (option1 easiest) :
Run "windows_launcher.bat" file as administrator. It will open a powershell window as well as a browser tab. If the browser tab doesn't show anything after 10sec, refresh it. 


### Run in powershell/terminal (option2) :
> docker run --rm -p 5288:5288 `
    -e SHINY_PORT=5288 `
    -e SHINY_ROOT_PATH=/browse `
    -e SHINY_ROOT_NAME=home `
    --mount type=bind,source="C:\Users\username\Documents",target=/browse `
    --mount type=bind,source="C:\shiny_out",target=/out `
    qtea1/bulktools

nb : You should edit the mount paths so that they correspond to your device's.


## FOR MACOS (May not work for everyone) :
### Click style (option1 easiest) :
Double-click "macos_launcher.command". It opens a Terminal window and a browser tab. If the browser tab is empty after 10sec, refresh it.

First launch : macOS may block the script ("unidentified developer"). Right-click the file > Open > Open (once), or run in a terminal :
> xattr -d com.apple.quarantine macos_launcher.command

If double-clicking does nothing, make it executable once :
> chmod +x macos_launcher.command




nb : To build the docker use this type of cmd (optionnal) :
> docker build -t qtea1/bulktools .
