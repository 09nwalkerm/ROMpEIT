CI/CD Pipelines
===============

Should you wish to adapt the repository, it is recommended that the user sets up a continuous integration (CI)/continuous development (CD) pipeline.
A suite of tests has been written in the :code:`tests/` folder to make sure the code is behaving as it should.
In this section we describe a Gitlab-based method for setting up your own CI/CD pipeline.

Gitlab
------

Once the repository is downloaded, link it to a repository on a Gitlab instance. This can be done using the following line, replacing the necessary usernames and addresses.

.. code-block:: bash

   git remote add gitlab git@gitlab.com:username/ROMpEIT.git

Now, add a .gitlab-ci.yml file to the root of the repository. An example is provided in the :code:`setup/CI/` folder.
This file triggers a CI/CD pipeline on Gitlab that pulls a MATLAB-based Docker image and runs the container.
In the container, a suite of tests from the :code:`tests/` folder are executed.

Before pushing this change to Gitlab however, the MATLAB license file variable :code:`MLM_LICENSE_FILE` needs to be changed to the address of the your instituations license server.
Otherwise the tests will fail.

The Docker image pulled is a custom one made with Statistics and Machine Learning Toolbox loaded with MATLAB. This image is :code:`09nwalkerm:matlab:v2`.
If you require more toolboxes in MATLAB then you must create your own publically available Docker container (on `Docker Hub <https://hub.docker.com/>`_) with the relevant toolboxes.
There is a Dockerfile in the :code:`setup/CI/` folder to get you started.
More info on setting up custom MATLAB containers can be found `here <https://github.com/mathworks-ref-arch/matlab-dockerfile>`_.


Github Mirroring
----------------

If you have forked the Github repository of this toolbox and wish to keep pushing changes to Github, then Gitlab can still be used to run tests automatically through setting up a pull mirror on you Gitlab instance.
We recommend following the instructions laid out in `this guide <https://docs.gitlab.com/ee/ci/ci_cd_for_external_repos/github_integration.html>`_ to set up this pipeline.  
A brief overiew is as follows: On your Gitlab project, go to Settings > Repositories > Mirroring Repositories, and delete any existing mirror, adding the HTTPS URL from the Github project but with the addition of 'YOURUSERNAME@' right before 'github.com'. Add a Personal Access Token from Github and check the first 2 boxes to allow CI events.
