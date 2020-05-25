# MITE
MITE (Matrix of Intensity by Time of Elution)-based, phyloproteomic analyses

Steps to install and run this project:

[Download and install docker](https://docs.docker.com/get-docker/) 

You might want to [add your user to the Docker group](https://docs.docker.com/engine/install/linux-postinstall/) in order to run the docker commands without sudo.

In the project root folder, run the following commands:
```console
foo@bar:~$ docker image build -t mite .
```

```console
foo@bar:~$ docker container run -it --name mite -v $(pwd):/user/src/mite mite
```

Now you should have a console inside the project container. Run the following command to run the project pipeline.
```console
root@bar:/user/src/mite# cd src/mite/ && python pipeline.py --nproc 4 --window_height 12 --window-width 12
```

To exit the container console, run the exit command:
```console
root@bar:/user/src/mite# exit
```

In order to remove the container and image from your system, run the following command:
```console
foo@bar:~$ docker container rm mite && docker image rm mite
```

