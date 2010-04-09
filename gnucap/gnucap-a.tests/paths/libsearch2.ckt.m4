* test for library search 

v1 n1 0 dc 1
.options includepath="PATH1:PATH2"
.include ./lib_file2.ckt
.include ./lib_file3.ckt

.options
.list
.print op v(n*)
.op
.end
