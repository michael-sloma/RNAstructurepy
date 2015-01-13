 %module RNAstructure
 %{
 /* Includes the header in the wrapper code */
 #include "../RNA_class/RNA.h"
 #include "../RNA_class/Dynalign_object.h"
 #include "../RNA_class/Multilign_object.h"
 #include "../RNA_class/Oligowalk_object.h"
 #include "../RNA_class/ProbScan.h"
 %}
 
 /* Parse the header file to generate wrappers */
 %include "../RNA_class/RNA.h"
 %include "../RNA_class/Dynalign_object.h"
 %include "../RNA_class/Multilign_object.h"
 %include "../RNA_class/Oligowalk_object.h"
 %include "../RNA_class/ProbScan.h"