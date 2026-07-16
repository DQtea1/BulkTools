docker build -t qtea1/bulktools .


# WITH DOCKER : 

## 1) Install Docker Desktop
Download docker desktop then open it :

https://docs.docker.com/desktop/setup/install/windows-install/

## 2) Pull the docker image from Dockerhub :
Run in the terminal :
> docker pull qtea1/bulktools:latest

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


## FOR MACOS (two options) :
### Click style (option1 easiest) :
Double-click "macos_launcher.command". It opens a Terminal window and a browser tab. If the browser tab is empty after 10sec, refresh it.

nb (first launch) : macOS may block the script ("unidentified developer"). Right-click the file > Open > Open (once), or run in a terminal :
> xattr -d com.apple.quarantine macos_launcher.command

nb : if double-clicking does nothing, make it executable once :
> chmod +x macos_launcher.command

### Apple Silicon (M1 / M2 / M3 / M4) :
The image is built for Intel/amd64 only, so on Apple Silicon Docker runs it through Rosetta emulation.
1. In Docker Desktop > Settings > General, enable "Use Rosetta for x86/amd64 emulation".
2. The launcher already passes "--platform linux/amd64". If you run docker by hand, add that flag too (see below).

### Run in terminal (option2) :
> docker run --rm -p 5288:5288 \
    --platform linux/amd64 \
    -e SHINY_PORT=5288 \
    -e SHINY_ROOT_PATH=/browse \
    -e SHINY_ROOT_NAME=home \
    --mount type=bind,source="$HOME",target=/browse \
    --mount type=bind,source="$HOME/shiny_out",target=/out \
    qtea1/bulktools
nb : "--platform linux/amd64" is required on Apple Silicon and harmless on Intel Macs. Edit the mount paths if you want other folders.


## FOR LINUX (doesnt work yet):
> docker run --rm -p 5288:5288 \
    -e SHINY_PORT=5288 \
    -e SHINY_ROOT_PATH=/browse \
    -e SHINY_ROOT_NAME=home \
    -v /home/user:/browse \
    -v /home:/out \
    qtea1/bulktools
