
# WITHOUT DOCKER :
1. Set yourself in the same directory as the "app.R" file.

2. Run in R terminal :
> shiny::runApp(host = "0.0.0.0", port = 5288, launch.browser = FALSE)

3. Open in browser :
> http://localhost:5288

4. Follow the tutorial in each module


# WITH DOCKER (recommanded) : 

## Install Docker Desktop
Download docker desktop then open it :

https://docs.docker.com/desktop/setup/install/windows-install/

## Pull the docker image from Dockerhub :
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


## FOR LINUX (doesnt work yet):
> docker run --rm -p 5288:5288 \
    -e SHINY_PORT=5288 \
    -e SHINY_ROOT_PATH=/browse \
    -e SHINY_ROOT_NAME=home \
    -v /home/user:/browse \
    -v /home:/out \
    qtea1/bulktools
