## GTF/GFF parser ##
Given a list of gene/attribute ids or names and a GTF or GFF formatted file, parse the input file to find entries with the names or ids in attribute section. The input list of ids/names must be one per line. Supports reading gzipped (.gz) files.

### Dependencies ###
* gcc (version >=4.8 for c++11 standard support)
* [boost](http://www.boost.org/) library ( for gzip support)
	* boost-dev
	* boost-iostreams-dev

#### Boost library installation ####
On a Debian/Ubuntu machine install boost libraries using:  
`sudo apt-get install libboost-dev libboost-iostreams-dev`  
For other distros and unix variants, see [boost library installation](http://www.boost.org/doc/libs/1_61_0/more/getting_started/unix-variants.html)

### Installation ###
Clone the repository to your local disk:  
`git clone https://github.com/sudeepsahadevan/gtf_parser.git`  
Change to directory and make  
`cd gtf_parser`  
`make`  
The executable `gtf_parser` will be created in the folder `bin`

### Usage ###
Call the executable as:  
`/path/to/gtfparser --ids /path/to/ids --gff /path/to/gtf > gff_filtered`  
For more usage options see:  
`/path/to/gtfparser -h` 
