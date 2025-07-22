CI/CD Pipelines
===============

Should you wish to adapt the repository, it is recommended that the user sets up a continuous integration (CI)/continuous development (CD) pipeline.
A suite of tests has been written in the :code:`tests/` folder to make sure the code is behaving as it should.
In this section we describe a GitLab-based method for setting up your own CI/CD pipeline.

GitLab
------

Once the repository is downloaded, link it to an empty repository on a GitLab instance. Presuming the necessary permissions are set up, this can be done using the following lines, replacing the necessary usernames and addresses.

.. code-block:: bash

   git remote add gitlab git@gitlab.com:username/ROMpEIT.git
   git push gitlab --all

Now, add a .gitlab-ci.yml file to the root of the repository. An example is provided in the :code:`setup/CI/` folder.
This file triggers a CI/CD pipeline on GitLab that pulls a MATLAB-based Docker image and runs the container.
In the container, a suite of tests from the :code:`tests/` folder are executed.

Before pushing this change to GitLab however, the MATLAB license file variable :code:`MLM_LICENSE_FILE` needs to be changed to the address of the your instituations license server.
Otherwise the tests will fail.

The Docker image pulled is a custom one made with Statistics and Machine Learning Toolbox loaded with MATLAB. This image is :code:`09nwalkerm:matlab:v2`.
If you require more toolboxes in MATLAB then you must create your own publically available Docker container (on `Docker Hub <https://hub.docker.com/>`_) with the relevant toolboxes.
There is a Dockerfile in the :code:`setup/CI/` folder to get you started.
More info on setting up custom MATLAB containers can be found `here <https://GitHub.com/mathworks-ref-arch/matlab-dockerfile>`_.


GitHub Mirroring
----------------

If you have forked the GitHub repository of this toolbox and wish to keep pushing changes to GitHub, then GitLab can still be used to run tests automatically through setting up a pull mirror on you Gitlab instance.
We recommend following the instructions laid out in `this guide <https://docs.gitlab.com/ee/ci/ci_cd_for_external_repos/GitHub_integration.html>`_ to set up this pipeline.  
