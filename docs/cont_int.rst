CI/CD Pipelines
===============

Should you wish to fork and adapt the repository, it is recommended that the user sets up a continuous integration (CI)/continuous development (CD) pipeline.
A suite of tests has been written in the :code:`tests/` folder to make sure the code is behaving as it should.
In this section we describe a Gitlab based method for setting up your own CI/CD pipelines depending on your available resources.

Gitlab
------

- Please follow `this guide <https://docs.gitlab.com/ee/ci/ci_cd_for_external_repos/github_integration.html>`_ to set up a repo with github integration and appropriate webhook on this repo.

- If you require more toolboxes in matlab then you must create your own publically available (from `Docker Hub <https://hub.docker.com/>`_) docker container with the relevant toolboxes.

- More info on setting up custom matlab containers can be found `here <https://github.com/mathworks-ref-arch/matlab-dockerfile>`_.
  
- There is a Dockerfile in the :code:`setup/CI/` folder to get you started. 

- It's possible GitLab sets up the Repository Mirroring incorrectly. On GitLab project, go to Settings > Repositories > Mirroring Repositories, and delete the existing mirror, adding the HTTPS URL from the github project but with the addition of 'YOURUSERNAME@' right before 'github.com'. Add Personal Access Token from GitHub and check the first 2 boxes to allow CI events.

