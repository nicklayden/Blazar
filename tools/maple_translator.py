# -*- coding: utf-8 -*-
"""
    Class made for translating maple CodeGeneration optimized functions into
    cpp functions for use with my root finding algorithm in cpp.
    
    author: Nick Layden
"""
import os

class CreateFunction(object):
    
    def __init__(self, file_,ret_type_="mp_type",params=["z","eta"], outfile_=None,idir_=""):
        # File names only
        self.file_ = file_
        self.outfile_ = outfile_
                
        # After this many variables defined, wrap to new line for definitions.
        self.decl_wrap_length = 20
        self.return_type = ret_type_
        self.func_parameters = params
        self.ext = ".cpp"
        self.odir_ = "./"
        self.idir_ = idir_
        # Construct the file objects from file names above
        self.input_file = self.OpenInputFile()
        self.output_file = self.OpenOutputFile()
    
        # Find all variables necessary
        self.var_list = self.FindVariables()
    
        self.ConstructFunction()

        # Closing all files after constructfunction has finished.
        self.input_file.close()
        self.output_file.close()
    
    def OpenInputFile(self):
        print("Opening file: {}".format(self.idir_ + self.file_))
        return open(self.idir_+self.file_)
    
    def OpenOutputFile(self):
        if self.outfile_ == None:
            outfile_name = self.odir_ + self.file_.split(".")[0] + self.ext
            print("Creating file: {}".format(outfile_name))
            return open(outfile_name,"w+")
        else:
            print("Creating file: {}".format(self.outfile_.split(".")[0] + self.ext))
            return open(self.odir_ + self.outfile_.split(".")[0] + self.ext, "w+")
        
    def CountLines(self):
        i=0
        for line in self.input_file:
            i+=1
        return print(i)
    
    def FindVariables(self):
        var_list = []
        for line in self.input_file:
            var = line.split("=")[0]
            var_list.append(var)
        return var_list
    
    def ConstructVariableList(self):
        '''
            Takes the maple output file of variables, creates all necessary
            variable declarations and puts them in the output file.
            
            The list comp makes sublists of the total variable list of size "skip"
            or the remaining number of elements if there are less than skip vars
            in the list. Then joins each sublist to write into the out 
            file as a line
            var_line is the line of variables truncated to size "skip"
        '''
        skip = self.decl_wrap_length
        var_layered = [self.var_list[i:i+skip] for i in range(0,len(self.var_list),skip)] 
        for i in range(len(var_layered)):
            var_line = "\t" + self.return_type + " " + ", ".join(var_layered[i]) + ";\n"
            self.output_file.write(var_line)

        
    def ConstructFunctionHeader(self):
        ''' 
            Create the cpp function from the maple output, give it same name
            as the input file name.
            
            Using python string fuckery to create the proper argument list for 
            the function inputs.
        '''
        FuncName = self.file_.split(".")[0]
        FuncArgsList = [self.return_type + " " + arg for arg in self.func_parameters]
        FuncArgsStr = ", ".join(FuncArgsList)
        FuncHeader = self.return_type + " " + FuncName + "(" + FuncArgsStr + ")\n"
        self.output_file.write(FuncHeader)
        
    def ConstructFunctionInitialBrace(self):
        self.output_file.write("{\n")
        
    def ConstructFunctionFinalBrace(self):
        self.output_file.write("\n}\n")
        
    def ConstructFunctionBody(self):
        self.output_file.write("\n\n")
        with open(self.idir_+self.file_) as f:
            for line in f:
                self.output_file.write("\t" + str(line))
            
    def ConstructReturnStatement(self):
        self.output_file.write("\n\treturn {};".format(self.var_list[-1]))
    
    def ConstructFunction(self):
        print("Building function...")
        print("Creating Header")
        self.ConstructFunctionHeader()
        print("Creating opening brace")
        self.ConstructFunctionInitialBrace()
        print("Creating variable declarations")
        self.ConstructVariableList()
        print("Creating function body")
        self.ConstructFunctionBody()
        print("Creating return statement")
        self.ConstructReturnStatement()
        print("Creating closing brace")
        self.ConstructFunctionFinalBrace()
            
            
            
            
# CreateFunction("app_horizon.txt", outfile_="app_horizon",idir_="../mfunctions/")

# mdir = "../mfunctions/"
# for item in os.listdir(mdir):
#     print(item)
#     print(str(item.split(".")[0]))
#     CreateFunction(item, outfile_= str(item.split(".")[0]),idir_=mdir)  
        
        
        
CreateFunction("eta_crunch.txt",outfile_="eta_crunch",idir_="../mfunctions/")        
CreateFunction("eta_bang.txt",outfile_="eta_bang",idir_="../mfunctions/")        
#CreateFunction("R_example1.txt",outfile_="R_example1",idir_="../mfunctions/")        

        
        
        
        
        
        
        
        
        
    



















