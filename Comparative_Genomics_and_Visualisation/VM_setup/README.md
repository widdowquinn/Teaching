# README.md - Virtual Machine Setup

For this course, we will use a VirtualBox CentOS 6.4 virtual image.

# 1. Obtaining and Installing VirtualBox

## a. Installing VirtualBox

VirtualBox is virtualisation software. It allows you, effectively, to run one computer operating system inside another. We will use it on this course to host a virtual machine (VM) running the CentOS 6.4 operating system. The VirtualBox site can be found at <https://www.virtualbox.org/>. 

![VirtualBox homepage](images/virtualbox1.png?raw=True =300x)

You should use the **Downloads** link on the left of the home page, and select the appropriate version of the software for your operating system.

![VirtualBox download page](images/virtualbox2.png?raw=True =300x)

Download the appropriate installation package, and follow the instructions to install the software on your machine.

## b. Installing a Virtual Machine (VM)

VirtualBox has the ability to run a very large range of guest operating systems. A sample of these, with download links can be seen at <http://virtualboxes.org/images/>.

![VirtualBoxes download page](images/virtualboxes1.png?raw=True =300x)

This course was written and developed for the CentOS 6.4 i386 target OS. A suitable VM, `CentOS-6.4-i386-Gnome.ova`, may be downloaded from <http://virtualboxes.org/images/centos/>.

Additionally, the course has been tested on an Ubuntu 12.10 VM. A suitable VM, `ubuntu-12.10-desktop-i386`, can be downloaded directly from <http://sourceforge.net/projects/virtualboximage/files/Ubuntu%20Linux/12.10/ubuntu-12.10-desktop-i386.7z>

### CentOS 6.4

Once the `CentOS-6.4-i386-Gnome.ova` VM has downloaded, double-clicking on it should start the VirtualBox application:

![CentOS import](images/centos1.png?raw=True =300x)

Click on `Import`, and the VM should be imported, and show you details:

![CentOS import result](images/centos2.png?raw=True =300x)

Initially, the VM's RAM will be set to 512MB. This is too small for some of the activities, so change it to something larger (I set 4GB), by clicking on `Settings -> System` and moving the slider:

![CentOS settings](images/centos3.png?raw=True =300x)

Then double-click the CentOS-6.4 image at the left of the VM manager to start the machine.


### Ubuntu 12.10

The `ubuntu-12.10-desktop-i386` VM downloads as the file `ubuntu-12.10-desktop-i386.7z`, which is in `7zip` format. You may need to download a suitable unarchiving utility from `http://www.7-zip.org/` or elsewhere to open this file.

On extracting the contents of the file, you should see a new directory has been created, called `ubuntu-12.10-desktop-i386`. You should move this directory to the location in which VirtualBox stores its VMs. (On my Mac laptop, this is in `~/VirtualBox\ VMs`).

You can then start the VM for the first time, by double-clicking on the file `ubuntu-12.10-desktop-i386.vbox`, in the `ubuntu-12.10-desktop-i386` directory. You should see that the new VM is present in VirtualBox:

![VirtualBox Ubuntu details](images/ubuntu1.png?raw=True =300x)

To start the VM, click on the `Start` arrow at the top of the window, or double-click on the link at the left.

**NOTE:** With the downloaded VM that I tried, the default keyboard layout was Italian. You will probably need to change this (there is a tutorial on how to do this at <http://www.youtube.com/watch?v=ByR3uTWe1yI>).


# 2. Installing course materials and bioinformatics software

The disk image does not by default contain all the software and tools required for the activities and exercises in this workshop, and the rest of this Markdown document describes the process of preparing the VM for use.

## a. The Easy Way

To begin, start a Terminal session in the VM. 

You will need first to install `git`, then download the course materials, and the `i-ADHoRe` software. Finally, you will need to run a `Makefile` to install bioinformatics software for the course.

### Install `git`

To download the course materials, you will need to install `git`. This is not present by default on either VM. To install `git`, follow the instructions below for your operating system.


#### CentOS
Become the root user, by issuing:

```
$ su -
<enter root password>
```

and then install `git` with:

```
$ yum -y install git
```

You will be asked to confirm this, and you should. Then use `Crtl-D` (i.e. press and hold the `Ctrl` key, and press `d`) to return to your original user shell - or start another terminal so you can work in parallel as root and the standard user.

#### Ubuntu

Issue the command:

```
$ sudo apt-get -y install git
```

When this is complete, you will be ready to 

### Download course materials

The instructions are the same for both CentOS and Ubuntu VMs.

Change to your home directory with the command:

```
$ cd
```

and issue the following command to download all the comparative genomics workshop teaching materials, and VM setup script, to your VM in the current directory:

```
$ git clone https://github.com/widdowquinn/Teaching.git
```

When this is complete, the course materials will have been downloaded to the subdirectory `Teaching`

### Download `i-ADHoRe`

Once the course materials have downloaded, use the Firefox web browser (either use the system menus, the icon, or issue `firefox` at the command line) to go to the website <http://bioinformatics.psb.ugent.be/software/details/i--ADHoRe>.

If you are using the **CentOS** VM, you will have to install Firefox first, by becoming root, and using `yum`:

```
$ su -
$ yum -y install firefox
```

and using `Ctrl-D` to return to being the normal `centos` user.

![i-ADHoRe website](images/iadhore1.png?raw=True =300x)

Click on the "Download Software" link, and agree to the terms and conditions of the agreement, for `i-ADHoRe 3.0`. When asked, choose to save the file (rather than open with the Archive Manager).

![i-ADHoRe agreement](images/iadhore2.png?raw=True =300x)

This will place the `i-ADHoRe-3.0.01.tar.gz` file in your `~/Downloads` directory. Move this to the appropriate location by issuing the command:

```
$ mv ~/Downloads/i-ADHoRe-3.0.01.tar.gz ~/Teaching/Comparative_Genomics_and_Visualisation/VM_setup
```

### Install course software and activity data

Now, change directory to the `VM_setup` directory:

```
$ cd ~/Teaching/Comparative_Genomics_and_Visualisation/VM_setup
```

The installation for activity data, and for each application and tool (where possible) is handled by a `Makefile`. There are two `Makefile`s: one for CentOS, and one for Ubuntu, as they use different package managers. Be sure to follow the instructions for your choice of operating system:

#### CentOS

Become root again with:

```
$ su -
<enter root password>
```

and copy the `Makefile_CentOS` file with:

```
cp Makefile_CentOS Makefile
```

Then issue the command `make all`:

```
$ make all
```

Then use `Crtl-D` to become the user `centos` again


#### Ubuntu

copy the `Makefile_Ubuntu` file with:

```
cp Makefile_Ubuntu Makefile
```

Then issue the command `make all` with root privileges:

```
$ sudo make all
```



## b. As individual packages

## Python

`Python` is available by default on CentOS 6. Although it is version 2.6, all activities and exercises should run normally. The Ubuntu VM provides Python 2.7.3 by default.

## GCC

`GCC` is a compiler, necessary for some of the subsequent download and installation steps outlined below. To install, as root, issue the following:

#### CentOS

```
$ yum -y install gcc gcc-c++
```

#### Ubuntu

```
$ apt-get -y install gcc g++
```


## FireFox

`FireFox` is a browser, and required to use the iPython notebook, employed in examples and activities. It is not installed on the CentOS 6.4 VM by default, but is present on Ubuntu 12.10. To install, issue as root:

#### CentOS

```
$ yum -y install firefox
```

## iPython and other Python libraries

`iPython` is an interactive environment for Python, and is required to run several of the activities. It has dependencies on other libraries. To install as root, issue :

#### CentOS

`python-pip` is in the `EPEL` repository, which is not installed by default on CentOS:

```
$ yum -y install wget
$ wget http://mirror-fpt-telecom.fpt.net/fedora/epel/6/i386/epel-release-6-8.noarch.rpm
$ rpm -ivh epel-release-6-8.noarch.rpm
$ yum -y install zeromq scipy sympy python-nose python-pip python-matplotlib
$ pip -y install ipython[all]
```

#### Ubuntu

```
$ apt-get -y install libzmq-dev python-scipy python-sympy python-nose python-pip python-matplotlib
$ pip -y install ipython[all]
```


## ReportLab

`ReportLab` is a graphics library for Python, and is required for Biopython visualisation modules. To install, issue the following as root:

#### CentOS

```
$ yum -y install python-reportlab python-imaging
```

#### Ubuntu

```
$ apt-get -y install python-reportlab python-imaging
```


## BioPython

`Biopython` is a Python bioinformatics library, required for several activities and examples. To install, issue the following as root:

#### CentOS/Ubuntu

```
$ pip install biopython
```

## Pandas

`Pandas` is a Python module for dataframe manipulation. It is required in some examples and activities. It needs a more recent version of `numpy` than provided by CentOS. To install, issue as root:

#### CentOS/Ubuntu

```
$ pip install numpy --upgrade
$ pip install pandas
```


## R and RStudio

`R` is a statistical programming environment, much used in bioinformatics. It is required for some activities and examples, as is the IDE `RStudio`. As root, issue:

#### CentOS

```
$ yum -y install R
$ wget http://download1.rstudio.org/rstudio-0.98.501-i686.rpm
$ yum install --nogpgcheck rstudio-0.98.501-i686.rpm
$ rm rstudio-0.98.501-i686.rpm
```

#### Ubuntu

```
$ add-apt-repository "deb http://cran.rstudio.com/bin/linux/ubuntu quantal/"
$ apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
$ apt-get update
$ apt-get -y --force-yes install r-base r-base-dev
$ wget http://download1.rstudio.org/rstudio-0.98.501-i386.deb
$ dpkg -i rstudio-0.98.501-i386.deb
$ rm rstudio-0.98.501-i386.deb
```

#### CentOS/Ubuntu

It will also be necessary to install further R libraries. You can start `RStudio` and use the GUI, or install from the `R` command-line. The required packages are:

* **ggplot2**
* **hash**
* **gridExtra**
* **reshape2**
* **shiny**

And to install **BioConductor** (to use **Biostrings**), we could install using:

```
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("Biostrings")
```

at the R command-line.

Alternatively, once `R` is installed, you could run the provided script:

```
$ Rscript install_R_packages.R
```

## Mercurial

In order to get a working copy of `RPy2` we have to get the source using Mercurial. To install, issue as root:

#### CentOS

```
$ yum -y install hg
```

#### Ubuntu

```
$ apt-get -y install mercurial meld
```

## RPy2

`RPy2` is a Python module for interacting with the `R` software. To install, issue the following as root to get the source from BitBucket and install version 2.1 (the most recent that works on the VM provided):

Then use **Mercurial** to fetch the source code, and install as root:

```
$ hg clone https://www.bitbucket.org/lgautier/rpy2
$ cd rpy2
$ hg update -C version_2.1.x
$ python setup.py build
$ python setup.py install
$ cd ..
$ rm -rf rpy2
```


## BLAST+ and Legacy BLAST

`BLAST+` is required for activities and exercises, and the legacy version of `BLAST` is required for `JSpecies` and some activities and exercises. To install, issue the following as root:

```
$ wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.2.29+-ia32-linux.tar.gz
$ tar -zxvf ncbi-blast-2.2.29+-ia32-linux.tar.gz
$ mv ncbi-blast-2.2.29+ /opt/
$ rm ncbi-blast-2.2.29+-ia32-linux.tar.gz
$ wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/2.2.26/blast-2.2.26-ia32-linux.tar.gz
$ tar -zxvf blast-2.2.26-ia32-linux.tar.gz
$ mv blast-2.2.26 /opt/
$ rm blast-2.2.26-ia32-linux.tar.gz
```

The paths `/opt/blast-2.2.26/bin` and `/opt/ncbi-blast-2.2.29+/bin` then need to be added to your `${PATH}`. This can be done using the `path_extensions.txt` file:

```
$ cat path_extensions.txt >> /home/centos/.bashrc
```

## MUMmer

`MUMmer` is a sequence alignment package, required for several activities and exercises. To install, issue the following as root:

```
$ wget http://downloads.sourceforge.net/project/mummer/mummer/3.23/MUMmer3.23.tar.gz
$ tar -zxvf MUMmer3.23.tar.gz
$ mv MUMmer3.23 /opt/
$ cd /opt/MUMmer3.23
$ make install
$ cd
$ rm MUMmer3.23.tar.gz
```

Note that MUMmer has to be built **in-place** in `/opt/MUMmer3.23`, or else it cannot find support files, once it has been moved.

The path `/opt/MUMmer3.23` then needs to be added to your `${PATH}`. This can be done using the `path_extensions.txt` file:

#### Ubuntu

```
$ cat path_extensions.txt >> /home/centos/.bashrc
```


#### Ubuntu

```
$ cat path_extensions.txt >> /home/ubuntu/.bashrc
```

## i-ADHoRe

`i-ADHoRe` is a tool for identifying syntenous and collinear genome features, required for activities. We cannot distribute or automatically download the source, so this will need to be downloaded separately from <http://bioinformatics.psb.ugent.be/software/details/i--ADHoRe>.

To install on Linux, issue:

```
tar -zxf i-adhore-3.0.01.gz
mkdir -p i-adhore-3.0.01/build
cd i-adhore-3.0.01/build
cmake .. -DCMAKE_INSTALL_PREFIX=/opt/i-adhore_3.0
make
make install
```