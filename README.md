
# WITHOUT DOCKER :
1. Move yourself in the same directory as the "app.R" file.

2. Run in R terminal :
> shiny::runApp(host = "0.0.0.0", port = 5288, launch.browser = FALSE)

3. Open in browser :
> http://localhost:5288



# WITH DOCKER (recommanded) : 

## Install Docker Desktop
https://docs.docker.com/desktop/setup/install/windows-install/

## Build the docker yourself (option1) :
Set yourself in the same directory as "Dockerfile" then run in the terminal :
> docker build -t bulktools .

## Pull the docker image from Dockerhub (option2 recommanded):
Run in the terminal :
> docker pull qtea1/bulktools:latest

## RUN FOR WINDOWS (two options) : 
### Click style (option1 easiest) :
Right click and execute as adminstrator on "windows_launcher.bat".

nb : Default root for selecting an output path will be "C:\Users" so you might want to edit it to a customed path.

### Run in powershell/terminal (option2) :
> docker run --rm -p 5288:5288 `
    -e SHINY_PORT=5288 `
    -e SHINY_ROOT_PATH=/browse `
    -e SHINY_ROOT_NAME=home `
    --mount type=bind,source="C:\Users\Quentin\Documents",target=/browse `
    --mount type=bind,source="C:\shiny_out",target=/out `
    qtea1/bulktools
nb : Default root for selecting an output path will be "C:\Users" so you might want to edit it to a customed path.


## FOR LINUX (doesnt work yet):
> docker run --rm -p 5288:5288 \
    -e SHINY_PORT=5288 \
    -e SHINY_ROOT_PATH=/browse \
    -e SHINY_ROOT_NAME=home \
    -v /home/quentin:/browse \
    -v /home:/out \
    qtea1/bulktools
