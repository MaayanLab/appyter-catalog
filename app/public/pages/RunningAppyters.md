# Running appyters

Appyters are served by us publicly, but can also be run by you on your local machine or in the cloud! There are several ways you can do this depending on what you're trying to run. The easiest way to get a fully working environment is to use the docker images we've already built for the appyters, these contain all the necessary dependencies you would need to run the appyter.

You can run the appyter application as a dockerized command line application which can allow you to
- serve a live jupyter notebook containing a notebook you've created *using* the public appyter platform to re-execute or modify
- serve it as a web-form as it appears on our system
- create or execute appyters as a command line application for instance as part of a computational workflow

## Step 1: Install docker
Download, install, and run Docker on your computer from [the official Docker website](https://www.docker.com/community-edition) if you do not yet have docker installed.

## Step 2. Collect necessary information to construct the appyter image
This information should be pre-populated if you came to this page from an appyter, alternatively you can enter it yourself.
