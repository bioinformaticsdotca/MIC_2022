---
layout: aws_page
permalink: /AWS_setup
title: AWS Setup
header1: Workshop Pages for Students
header2: AWS Setup
image: /site_images/Bioinfo_Logo.jpg
home: https://bioinformaticsdotca.github.io/high-throughput_biology_2017
---

### Logging into the Amazon Cloud During the Workshop  

* These instructions will **ONLY** be relevant in class, as the Cloud will not be accessible from home in advance of the class.
 
* On the cloud, we're going to use the default username: **ubuntu**
 

#### Logging in with ssh (Mac/Linux) <a id="maclinux"></a>
<p>
<h5> Logging in </h5>
</p>
* Make sure the permissions on your certificate are secure. Use chmod on your downloaded certificate:

```bash
 chmod 600 MIC.pem
```

* To log in to the node, use the -i command line argument to specify your certificate:

```bash
 ssh -i MIC.pem ubuntu@##.uhn-hpc.ca
```

(where ## is your assigned student number.)

* Everything created in your workspace on the cloud is also available by a web server on your cloud instance.  Simply go to the following in your browser:

```
 http://##.uhn-hpc.ca/
```

#### Logging in with Putty (Windows) <a id="windows"></a>
<p>
<h5> Logging in </h5>   
</p>
To configure Putty, start Putty and do the following:

* Fill in the "Host name" field with ##.uhn-hpc.ca (where ## is your assigned student number.)
 
<img src="https://github.com/bioinformaticsdotca/AWS_stuff/blob/master/Putty_Basic_Options.png?raw=true" alt="Basic Putty Options" class="center">

* In the left hand categories,under the Connection category choose Data.  In the auto-login username field write ***ubuntu***.

<img src="https://github.com/bioinformaticsdotca/AWS_stuff/blob/master/Putty_Data_Options.png?raw=true" alt="Putty Data Options" class="center"> 

* In the left hand categories, in the Connection category next to SSH click on the **+**. Click on Auth. In the private-key file for authentication field, hit browse and find the CBWHT.ppk certificate that you downloaded above.

<img src="https://github.com/bioinformaticsdotca/AWS_stuff/blob/master/Putty_Auth_Options.png?raw=true" alt="Putty Auth Options" class="center">

* In the left hand categories, click on Session.  In the Saved Sessions field write **Amazon node** and click save.

**Now that Putty is configured**, all you have to do is start putty and double-click on "Amazon node" to login.


#### File System Layout <a id="filesystem"></a>
<br>
When you log in, you'll notice that you have two directories: **CourseData** and **workspace**.

* The **CourseData** directory will contain the files that you'll need to complete your lab assignments.

* The **workspace** directory is where we will work on your labs. Files here will be accessible through your browser.

#### RStudio 
* RStudio server is installed on your instance. It is accessible in your browser
```
 http://##.uhn-hpc.ca:8080
```
(where ## is your assigned student number.) The username is ***ubuntu***. We will give you the password in Slack.
