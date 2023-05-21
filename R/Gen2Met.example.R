#Script: Gen2Met.example.R
#License: GPLv3 or later
#Modification date: 2023-05-03
#Written by: Valentino Palombo
#Contact: valentino.palombo@gmail.com
#Description: Copy the example files in the wd

Gen2Met.example<-function(type=NULL){
  #Copy files to current directory
  if (!(is.null(type))){
    if (type=="validation"){
      file.copy(system.file("data","Class1_annotated.txt", package = "Gen2Met"), "Class1_annotated.txt")
      file.copy(system.file("data","Class4_annotated.txt", package = "Gen2Met"), "Class4_annotated.txt")
    }else{
      file.copy(system.file("data","data.txt", package = "Gen2Met"), "data.txt")
    }
  }else{
    file.copy(system.file("data","data.txt", package = "Gen2Met"), "data.txt")
  }
}
