# ROMEG

## Continuous Integration

A breif guide for setting up CI with this repository for contributors and developers.

### Cardiff GitLab Instance

- Please follow [this guide](https://docs.gitlab.com/ee/ci/ci_cd_for_external_repos/github_integration.html) to set up a repo with github integration and appropriate webhook on this repo.
	- The webhook should start with "https://git.cardiff.ac.uk/api/v4/projects..." as this is a Cardiff University GitLab instance and not GitLab.com.

- Please create a new branch and edit the .gitlab-ci.yml file to run tests for only that branch.

- If you require more toolboxes in matlab then you must create your own publically available (from [DockerHub](https://hub.docker.com/)) docker container with the relevant toolboxes.
	- More info on setting up custom matlab containers can be found [here](https://github.com/mathworks-ref-arch/matlab-dockerfile)
	- There is a Dockerfile in the docker_context/ folder to get you started. 

- It's possible GitLab sets up the Repository Mirroring incorrectly. On GitLab project, go to Settings > Repositories > Mirroring Repositories, and delete the existing mirror, adding the HTTPS URL from the github project 
but with the addition of 'YOURUSERNAME@' right before 'github.com'. Add Personal Access Token from GitHub and check the first 2 boxes to allow CI events.

### Local Jenkins Server

- Download jenkins onto local machine and check jenkins is running: `sudo systemctl status jenkins`

- Jenkins should have created a user on Linux called 'jenkins'. Give this user a password, go to that userspace home directory and set up ssh config for allowing jenkins access to github so it can pull repo.

- Go to jenkins address on web browser (usually localhost:8080) and set up an account that isn't root. Create a Personal Access Token for this user and save it as you will need it for the github webhook.

- Set up a new project in jenkins and configure it to test the branch you have created.

- Download ngrok and start a webhook verified tunnel to expose the jenkins server port.

- Configure the setup/ci_setup file as necessary for whichever port you have as default for jenkins.

- Go to GitHub repo and create webhook with the address ngrok provided and the personal access token you made for the jenkins user.

- Create Deploy keys for github repo and save to local repo in relavent file (check this file matches .gitignore)
	- Now on each boot up of wsl you should be able to just run the ci_setup script and the webhook should be automatically updated.

### Matlab License

- For GitLab, the license should be left as is. The docker container is pointed to the Cardiff University Matlab License Server. This is however limited to version r2021a of Matlab.

- For Jenkins, the .bash_profile script of the jenkins user on Windows Subsystem for Linux needs to be editted to set a constant MAC address (eth0) that the Matlab license is linked to. You then need to go to 
the [Mathsworks License Center](https://uk.mathworks.com/?s_tid=gn_logo) and enter the static MAC address you used, activate the license and then download and save on the license path for Matlab in the jenkins userspace
of WSL.
