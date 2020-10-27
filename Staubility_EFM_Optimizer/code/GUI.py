from PyQt5.QtWidgets import *
import sys
from PyQt5 import QtGui,QtCore, QtWidgets
from PyQt5.QtCore import QRect
from functools import partial
from EFM import EFM_New

import time
import os


class PushButton:
    """
           A class used to implement Push button

           ...

           Attributes
           ----------
           button : QPushButton
               the actual push button from PyQt5 library

           Methods
           -------
           init(self)
               init the push button
           """
    def __init__(self,parent, name, location, Tooltip,function,arg1=None,arg2=None,arg3=None,arg4=None,arg5=None):
        """
        Parameters
        ----------
        parent : QGroupBox
            pointer to parent groupbox
        name : str
            The name of the button - will be displayed on the button
        location : 4dim location
            location of the button. suppose to get the next input:  QRect(10, 150, 120, 20)
        Tooltip : str
            description of the pushbutton functionality
        function : function
            link to the function that the pushbutton call to
        registers : list of register
            global registers list
        arg1,arg2,arg3,arg4,arg5:
            optiaonal input arguments to the callback function
        """
        self.button = QPushButton(parent=parent,text=name)
        self.button.setGeometry(location)
        self.button.setToolTip(Tooltip)
        if arg1==None:
            self.button.clicked.connect(partial(function))
        elif arg2==None:
            self.button.clicked.connect(partial(function, arg1))
        elif arg3==None:
            self.button.clicked.connect(partial(function, arg1,arg2))
        elif arg4 == None:
            self.button.clicked.connect(partial(function, arg1, arg2,arg3))
        elif arg5==None:
            self.button.clicked.connect(partial(function, arg1, arg2, arg3,arg4))
        else:
            self.button.clicked.connect(partial(function, arg1, arg2, arg3, arg4,arg5))

class FreeCheckButton:
    """
               A class used to implement Free Check button

               ...

               Attributes
               ----------
               button : QCheckButton
                   the actual check button from PyQt5 library

               Methods
               -------
               init(self)
                   init the check button
               """
    def __init__(self,parent, name, location, Tooltip,function=None,arg1=None,arg2=None,arg3=None,arg4=None,arg5=None):
        """
               Parameters
               ----------
               parent : QGroupBox
                   pointer to parent groupbox
               name : str
                   The name of the button - will be displayed near the button
               location : 4dim location
                   location of the button. suppose to get the next input:  QRect(10, 150, 120, 20)
               Tooltip : str
                   description of the checkbutton functionality
               function : function
                    link to the function that the pushbutton call to
               registers : list of register
                    global registers list
               arg1,arg2,arg3,arg4,arg5:
                    optiaonal input arguments to the callback function
               """
        self.button = QCheckBox(name, parent)
        self.button.setGeometry(location)
        self.button.setToolTip(Tooltip)
        self.button.setChecked(0)
        if function==None:
            return
        if arg1==None:
            self.button.clicked.connect(partial(function))
        elif arg2==None:
            self.button.clicked.connect(partial(function, arg1))
        elif arg3==None:
            self.button.clicked.connect(partial(function, arg1,arg2))
        elif arg4 == None:
            self.button.clicked.connect(partial(function, arg1, arg2,arg3))
        elif arg5==None:
            self.button.clicked.connect(partial(function, arg1, arg2, arg3,arg4))
        else:
            self.button.clicked.connect(partial(function, arg1, arg2, arg3, arg4,arg5))

class FreeComboBox:
    """
      A class used to implement ComboBox -using for decimal or str values
      this comboBox will NOT change values in the registers and fields structure
      it will change other values
      ...

      Attributes
      ----------
      button : QComboBox
          the actual ComboBox from PyQt5 library

      Methods
      -------
      init(self)
          init the ComboBox
      """
    def __init__(self, parent, items, location, Tooltip,function,value=None):
        """
           Parameters
           ----------
           parent : QGroupBox
               pointer to parent groupbox
           name : str
               The name of the button - not in use
           items: list
                list of strings to display
           location : 4dim location
               location of the button. suppose to get the next input:  QRect(10, 150, 120, 20)
           Tooltip : str
               description of the ComboBox functionality
           address : link to variable
               name of the variable to get the chosen value
           function : link to function
               callback function of this button
           """
        self.button = QComboBox(parent=parent)
        self.button.setGeometry(location)
        self.button.setToolTip(Tooltip)
        self.button.addItems(items)
        if value!=None:
            self.button.setCurrentText(value)
        self.button.currentIndexChanged.connect(partial(function, self.button))

class FreeLineEdit:
    """
              A class used to implement LineEdit
              this LineEdit will NOT change values in the registers and fields structure
              it will change other values
              ...

              Attributes
              ----------
              button : QLineEdit
                  the actual LineEdit from PyQt5 library

              Methods
              -------
              init(self)
                  init the LineEdit
              """
    def __init__(self, parent, location, Tooltip,value=None,function =None):
        """
               Parameters
               ----------
               parent : QGroupBox
                   pointer to parent groupbox
               location : 4dim location
                   location of the button. suppose to get the next input:  QRect(10, 150, 120, 20)
               Tooltip : str
                   description of the ComboBox functionality
        """
        self.button = QLineEdit(parent=parent)
        self.button.setGeometry(location)
        self.button.setToolTip(Tooltip)
        if value==None:
            self.button.setText("16")
        else:
            self.button.setText(value)
        if function!= None:
            self.button.editingFinished.connect(partial(function, self.button))






class App(QMainWindow):
    def __init__(self):
        super().__init__()
        self.title = "Staubility Tool"
        self.left = 100
        self.top = 100
        self.width = 800
        self.height = 600
        self.icon = QtGui.QIcon("In\\STAUBILITY2.png")
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        self.setWindowIcon(self.icon)

        # # Add tabs to widget
        self.table_widget = InitApp(self)
        self.setCentralWidget(self.table_widget)

        self.show()


class InitApp(QWidget):
    def __init__(self, parent):
        super(QWidget, self).__init__(parent)

        self.organism_name_list =["not_specified","e_coli","b_subtilis", "c_elegans","d_melanogaster",  "g_gallus", "h_sapiens", "m_musculus", "m_musculus_domesticus", "s_cerevisiae"]
        self.codon_methods = ["use_best_codon", "match_codon_usage", "harmonize_rca"]
        #
        self.input_path = "None"
        self.output_path = "None"
        self.consider_sites = False
        self.rec = "recA+"
        self.methylation_sites_path = r"In\topEnriched.313.meme.txt"
        self.number_of_sites = "All"
        self.optimize = False
        self.miniGC = 0.3
        self.maxiGC = 0.7
        self.method = self.codon_methods[0]
        self.organism_name = self.organism_name_list[0]
        self.indexes = None
        self.max_to_optimize = 10
        self.system_first_up = True


        self.main = parent
        self.layout = QVBoxLayout(self)

        # Initialize tab screen
        self.tabs = QTabWidget()
        self.tabs.resize(300, 200)

        # Initialize tab screen
        self.gsa = QTabWidget()


        # Add tabs area
        self.tabs.addTab(self.gsa, "Genomic Stability Analyzer")







        self.init_tab()


        # set layout
        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)

    def browse_path(self,input=True):
        if input:
            # load tx registes file
            options = QFileDialog.Options()
            #options |= QFileDialog.DontUseNativeDialog
            fileName = QFileDialog.getExistingDirectory(self, "Select Folder Of Fasta Files", "Folder",
                                                      options=options)
            if fileName:
                self.input_path = fileName
                self.all_gsa.input_path_line.button.setText(fileName)
        else:
            # load tx registes file
            options = QFileDialog.Options()
            #options |= QFileDialog.DontUseNativeDialog
            fileName = QFileDialog.getExistingDirectory(self, "Select Folder For Output File", "Folder",
                                                      options=options)
            if fileName:
                self.output_path = fileName
                self.all_gsa.output_path_line.button.setText(fileName)

    def change_sites_check(self):
        if self.all_gsa.sites.button.isChecked():
            QMessageBox.about(self, "User Warning", "This process will take a little longer, please use only if necessary")
        self.consider_sites = self.all_gsa.sites.button.isChecked()

    def change_rec_combo(self,button):
        self.rec = button.currentText()

    def change_sites_number_combo(self,button):
        number_of_sites = button.currentText()
        if number_of_sites=="ALL":
            self.number_of_sites = number_of_sites
        else:
            self.number_of_sites = int(number_of_sites)


    def change_organism_name(self,button):
        self.organism_name = button.currentText()


    def change_codon_method(self,button):
        self.method = button.currentText()


    def optomize_button(self):
        if self.all_gsa.optimize.button.isChecked():
            self.all_gsa.optimize.button.setChecked(False)
            if self.input_path == "None":
                QMessageBox.about(self, "User Warning",
                                  "Please Insert Input Path and press again on optimization button.")
            else:
                self.prop_indexes = EFM_New.advenced_function(input_folder=self.input_path)

                if len(self.prop_indexes) > self.max_to_optimize:
                    QMessageBox.about(self, "User Warning",
                                      "Oh No!\nThis Beta version support optimization of {} sequences or less.".format(self.max_to_optimize))
                elif not self.prop_indexes:
                    QMessageBox.about(self, "User Warning",
                                      "The system didn't find any files in the specified path\nPlease choose another path.")
                else:
                    self.dialog = QtWidgets.QDialog(self)
                    self.dialog.setWindowTitle("Optimization Settings")
                    self.dialog.setGeometry(100, 100, 750, 750)
                    #
                    self.dialog.optimize_group = QGroupBox("Optimization Parameters", self.dialog)
                    self.dialog.optimize_group.setGeometry(QRect(10, 10, 700, 100))
                    #
                    self.dialog.optimize_group.label_organ = QLabel("Organism Name:", self.dialog.optimize_group)
                    self.dialog.optimize_group.label_organ.setFont(QtGui.QFont('Arial', 12))
                    self.dialog.optimize_group.label_organ.setGeometry(QRect(10, 20, 150, 20))
                    self.dialog.optimize_group.organ_name = FreeComboBox(self.dialog.optimize_group, self.organism_name_list,
                                                                   QRect(215, 20, 120, 20), "Choose Organism Name",
                                                                   self.change_organism_name,self.organism_name)
                    #
                    self.dialog.optimize_group.label_codon = QLabel("Codon Optimization Method:",
                                                                     self.dialog.optimize_group)
                    self.dialog.optimize_group.label_codon.setFont(QtGui.QFont('Arial', 12))
                    self.dialog.optimize_group.label_codon.setGeometry(QRect(10, 50, 200, 20))
                    self.dialog.optimize_group.codon_method = FreeComboBox(self.dialog.optimize_group, self.codon_methods,
                                                             QRect(215, 50, 120, 20), "Choose Codon Optimizaion Method",
                                                             self.change_codon_method, self.method)
                    #
                    self.dialog.optimize_group.label_gc_content = QLabel("GC Content:", self.dialog.optimize_group)
                    self.dialog.optimize_group.label_gc_content.setFont(QtGui.QFont('Arial', 12))
                    self.dialog.optimize_group.label_gc_content.setGeometry(QRect(370, 20, 150, 20))
                    #
                    self.dialog.optimize_group.label_gc_content_min = QLabel("Min:", self.dialog.optimize_group)
                    self.dialog.optimize_group.label_gc_content_min.setFont(QtGui.QFont('Arial', 12))
                    self.dialog.optimize_group.label_gc_content_min.setGeometry(QRect(370, 50, 100, 20))
                    self.dialog.gc_min = FreeLineEdit(self.dialog.optimize_group, QRect(420, 50, 100, 20), "",
                                                       str(self.miniGC), self.change_gc_min)
                    #
                    self.dialog.optimize_group.label_gc_content_max = QLabel("Max:", self.dialog.optimize_group)
                    self.dialog.optimize_group.label_gc_content_max.setFont(QtGui.QFont('Arial', 12))
                    self.dialog.optimize_group.label_gc_content_max.setGeometry(QRect(540, 50, 100, 20))
                    self.dialog.gc_max = FreeLineEdit(self.dialog.optimize_group, QRect(590, 50, 100, 20), "",
                                                       str(self.maxiGC), self.change_gc_max)
                    #
                    self.dialog.optimize_group.label_files = QLabel("ORF Coding Regions",
                                                                     self.dialog)
                    self.dialog.optimize_group.label_files.setFont(QtGui.QFont('Arial', 13))
                    self.dialog.optimize_group.label_files.setGeometry(QRect(300, 110, 200, 30))
                    self.dialog.names = QLabel("File Name:", self.dialog)
                    self.dialog.names.setFont(QtGui.QFont('Arial', 12))
                    self.dialog.names.setGeometry(QRect(10, 150, 100, 20))
                    #
                    self.dialog.Loc = QLabel("Seq Num:", self.dialog)
                    self.dialog.Loc.setFont(QtGui.QFont('Arial', 12))
                    self.dialog.Loc.setGeometry(QRect(220, 150, 100, 20))
                    #
                    self.dialog.Start = QLabel("Start:", self.dialog)
                    self.dialog.Start.setFont(QtGui.QFont('Arial', 12))
                    self.dialog.Start.setGeometry(QRect(320, 150, 70, 20))
                    #
                    self.dialog.Stop = QLabel("Stop:", self.dialog)
                    self.dialog.Stop.setFont(QtGui.QFont('Arial', 12))
                    self.dialog.Stop.setGeometry(QRect(400, 150, 70, 20))
                    #
                    self.dialog.finish = PushButton(self.dialog, "Finished", QRect(370, 650, 100, 50),
                                                    "Save the data", self.close_and_save, True)
                    self.dialog.finish.button.setFont(QtGui.QFont('Arial', 12))
                    self.dialog.finish.button.setStyleSheet("background-color: green; font-weight: bold")
                    #
                    self.dialog.massege = QLabel(
                        "# Your inserted location must result in a sequence that is divisible by 3 and match the reading frame\n"\
                        "(target amino acid translation)." \
                        "The index convention: 1 is for start, X for end. For example, ATGCTG is (1,6).",
                        self.dialog)
                    self.dialog.massege.setFont(QtGui.QFont('Arial', 10))
                    self.dialog.massege.setGeometry(QRect(10, 550, 690, 80))
                    self.dialog.massege = QLabel(
                        "# You must press Finished button when you are done, to cancel please press X",
                        self.dialog)
                    self.dialog.massege.setFont(QtGui.QFont('Arial', 10,weight=QtGui.QFont.Bold))
                    self.dialog.massege.setGeometry(QRect(10, 580, 690, 80))

                    if self.system_first_up:
                        for i, obj in enumerate(self.prop_indexes):
                            ## Name
                            self.dialog.label = QLabel(self.prop_indexes[i][0], self.dialog)
                            self.dialog.label.setFont(QtGui.QFont('Arial', 8))
                            self.dialog.label.setGeometry(QRect(10, 190 + i * 30, 210, 20))
                            ## Index
                            self.dialog.label_index = QLabel(self.prop_indexes[i][1], self.dialog)
                            self.dialog.label_index.setFont(QtGui.QFont('Arial', 8))
                            self.dialog.label_index.setGeometry(QRect(250, 190 + i * 30, 50, 20))
                            ##
                            exec(
                                "self.dialog.start_{} = FreeLineEdit(self.dialog,QRect(320,190+i*30, 60, 20),'',str(1))".format(
                                    i))
                            exec(
                                "self.dialog.stop_{} = FreeLineEdit(self.dialog,QRect(400,190+i*30, 60, 20),'',str(self.prop_indexes[i][2]))".format(
                                    i))

                        self.close_and_save(close=False)
                        self.system_first_up = False
                    else:
                        i = 0
                        for obj in self.indexes:
                            ## Name
                            self.dialog.label = QLabel(obj[0], self.dialog)
                            self.dialog.label.setFont(QtGui.QFont('Arial', 8))
                            self.dialog.label.setGeometry(QRect(10, 190 + i * 30, 210, 20))
                            ## Index
                            self.dialog.label_index = QLabel(obj[1], self.dialog)
                            self.dialog.label_index.setFont(QtGui.QFont('Arial', 8))
                            self.dialog.label_index.setGeometry(QRect(250, 190 + i * 30, 50, 20))
                            ##
                            exec(
                                "self.dialog.start_{} = FreeLineEdit(self.dialog,QRect(320,190+i*30, 60, 20),'',str(self.indexes[obj][0]))".format(
                                    i))
                            exec(
                                "self.dialog.stop_{} = FreeLineEdit(self.dialog,QRect(400,190+i*30, 60, 20),'',str(self.indexes[obj][1]))".format(
                                    i))
                            i+=1


                    self.dialog.exec_()
        else:
            self.all_gsa.advanced.button.hide()

        self.optimize = self.all_gsa.optimize.button.isChecked()


    def change_gc_min(self,button):
        try:
            value = float(button.text())
            self.miniGC = value
        except:
            QMessageBox.about(self, "User Warning",
                              "GC Min must be a number")
            button.setText(str(self.miniGC))

    def change_gc_max(self,button):
        try:
            value = float(button.text())
            self.maxiGC = value
        except:
            QMessageBox.about(self, "User Warning",
                              "GC Max must be a number")
            button.setText(str(self.maxiGC))


    def open_user_guide(self):
        try:
            os.system("start EFM_User_Guide.pdf")
        except:
            QMessageBox.about(self, "User Warning",
                              "Oh No!\nCan't Find The UserGuide file. Please locate it in the main folder and keep the original name (UserGuide.txt)")

    def close_and_save(self, close=True):
        self.indexes = {}
        for i,obj in enumerate(self.prop_indexes):
            new_seq = []
            exec("new_seq.append(self.dialog.start_{}.button.text())".format(i)) # start
            exec("new_seq.append(self.dialog.stop_{}.button.text())".format(i))  # stop
            try:
                self.indexes[(obj[0],obj[1])] = [int(new_seq[0]),int(new_seq[1])]
            except:
                QMessageBox.about(self, "User Warning",
                                  "Error:\nFile: {} | Seq: {} - Indexes must be int.\n".format(obj[0],obj[1]))
                close = False

        if close:
            ret = self.check_optimization_parameters()
            if ret !="":
                QMessageBox.about(self, "User Warning",
                                  "There are some errors, please fix them and try again:\n\n"+ret)
            else:
                self.all_gsa.optimize.button.setChecked(True)
                self.all_gsa.advanced.button.show()
                self.dialog.close()



    def check_optimization_parameters(self):
        ret = ""
        ind = 0
        if self.miniGC<0 or self.miniGC>1:
            ret += "-Min GC must be between 0 to 1.\n"
        if self.maxiGC<0 or self.maxiGC>1:
            ret += "-Max GC must be between 0 to 1.\n"
        if self.miniGC > self.maxiGC :
            ret += "-Max GC must be grater than Min GC.\n"
        for key in self.indexes:
            if self.indexes[key][0]<1:
                ret+="-File: {} | Seq: {} - Start index must be grater than 0.\n".format(key[0],key[1])
            if self.indexes[key][1]> self.prop_indexes[ind][2]:
                ret += "-File: {} | Seq: {} - Stop index exceed the sequence length.\n".format(key[0], key[1])
            if self.indexes[key][0] > self.indexes[key][1]:
                ret += "-File: {} | Seq: {} - Start index must be grater the Stop index.\n".format(key[0], key[1])
            if (self.indexes[key][1] - self.indexes[key][0]+1)%3 !=0:
                ret += "-File: {} | Seq: {} - Sequence must be divided by 3.\n".format(key[0], key[1])
            ind+=1
        return ret



    def write_to_log(self,text):
        log = self.all_gsa.log.button
        log.setText(text)
        time.sleep(0.5)
        QtCore.QCoreApplication.processEvents()






    def start(self):
        if self.input_path == "None":
            QMessageBox.about(self, "User Warning",
                              "Error: Invalid Input Path!")
        elif self.output_path == "None":
            QMessageBox.about(self, "User Warning",
                              "Error: Invalid Output Path!")
        else:
            try:
                self.write_to_log("Starting Now - Please Wait Until The Process Will End!")

                ret = EFM_New.main(input_folder=self.input_path, output_path=self.output_path, compute_methylation=self.consider_sites,
                     num_sites=self.number_of_sites,
                     methylation_sites_path= self.methylation_sites_path, test=False, optimize=self.optimize,
                     miniGC=self.miniGC,
                     maxiGC=self.maxiGC, method=self.method, organism_name=self.organism_name, indexes=self.indexes)

                QMessageBox.about(self, "Message",ret)
                log = "Finished! you can find the results in the output path" if ret =="Success!" else ret
                self.write_to_log(log)
            except Exception as e:
                QMessageBox.about(self, "ERROR",
                                  "Oh No!\nSomething Went Wrong, Please Try Again!\n\nThe Error is:\n{}".format(e))



    def init_tab(self):
        self.all_gsa = QGroupBox("", self.gsa)
        self.all_gsa.setGeometry(QRect(0, 0, self.main.width, self.main.height))
        self.all_gsa.label = QLabel("Genomic Stability Analyzer",self.all_gsa)
        self.all_gsa.label.setFont(QtGui.QFont('Arial', 20))
        self.all_gsa.label.setGeometry(QRect(self.main.width/2-200, 20, 400, 35))
        #
        self.all_gsa.label1 = QLabel("Based On EFM Calculator",self.all_gsa)
        self.all_gsa.label1.setFont(QtGui.QFont('Arial', 15))
        self.all_gsa.label1.setGeometry(QRect(self.main.width/2-150, 45, 400, 35))
        #
        self.centralwidget = QtWidgets.QWidget(self.all_gsa)
        self.photo = QtWidgets.QLabel(self.centralwidget)
        self.photo.setGeometry(QtCore.QRect(600, 00, 170,170))
        self.photo.setPixmap(QtGui.QPixmap("In\\STAUBILITY2.png"))
        self.photo.setStyleSheet("border: 0")
        self.photo.setScaledContents(True)
        self.photo.setObjectName("photo")
        #

        self.all_gsa.label_input_path = QLabel("Input Sequence Path:", self.all_gsa)
        self.all_gsa.label_input_path.setFont(QtGui.QFont('Arial', 12))
        self.all_gsa.label_input_path.setGeometry(QRect(30, 100, 150, 20))
        self.all_gsa.input_path_line = FreeLineEdit(self.all_gsa,QRect(180, 100, 350, 20), "Insert Path To Fasta Files Folder","None")
        self.all_gsa.browse_input = PushButton(self.all_gsa,"Browse",QRect(550, 100, 50, 20),"Browse Input Path",self.browse_path,True)
        #
        self.all_gsa.label_output_path = QLabel("Output Path:", self.all_gsa)
        self.all_gsa.label_output_path.setFont(QtGui.QFont('Arial', 12))
        self.all_gsa.label_output_path.setGeometry(QRect(30, 140,150, 20))
        self.all_gsa.output_path_line = FreeLineEdit(self.all_gsa, QRect(180, 140, 350, 20), "Insert Path To Desired Output Folder","None")
        self.all_gsa.browse_output = PushButton(self.all_gsa,"Browse",QRect(550, 140, 50, 20),"Browse Output Path",self.browse_path,False)
        #
        self.all_gsa.label_sites = QLabel("Consider Methylation Sites:", self.all_gsa)
        self.all_gsa.label_sites.setFont(QtGui.QFont('Arial', 12))
        self.all_gsa.label_sites.setGeometry(QRect(30, 220,200, 20))
        self.all_gsa.sites = FreeCheckButton(self.all_gsa,"",QRect(240, 220,200, 20),"Organisms who express methylation in their genome are Saccharomyces cerevisiae… check with biologist",self.change_sites_check)
        #
        # self.all_gsa.label_rec =QLabel("recA:", self.all_gsa)
        # self.all_gsa.label_rec.setFont(QtGui.QFont('Arial', 12))
        # self.all_gsa.label_rec.setGeometry(QRect(30, 260,200, 20))
        # self.all_gsa.rec_combo = FreeComboBox(self.all_gsa,["recA+","recA-"],QRect(240, 260,100, 20),"Choose rec Type: Organisms with recA+ are… check with biologist",self.change_rec_combo)
        #
        self.all_gsa.label_number_sites = QLabel("Number Of Output Sites:", self.all_gsa)
        self.all_gsa.label_number_sites.setFont(QtGui.QFont('Arial', 12))
        self.all_gsa.label_number_sites.setGeometry(QRect(30, 260,200, 20))#QRect(30, 300, 200, 20))
        self.all_gsa.number_sites_combo = FreeComboBox(self.all_gsa,["ALL","5","10","20","50","100","150","200"],QRect(240, 260,100, 20),"Choose Number Of Output Sites",self.change_sites_number_combo)
        #
        self.all_gsa.label_optimize = QLabel("Optimize:", self.all_gsa)
        self.all_gsa.label_optimize.setFont(QtGui.QFont('Arial', 12))
        self.all_gsa.label_optimize.setGeometry(QRect(30, 300,200, 20))#QRect(30, 300, 200, 20))
        self.all_gsa.optimize = FreeCheckButton(self.all_gsa,"",QRect(240, 300,100, 20),"Optimize",self.optomize_button)
        #
        #### Optimize buttons sequence
        # self.all_gsa.optimize_group = QGroupBox("Optimization Parameters", self.all_gsa)
        # self.all_gsa.optimize_group.setGeometry(QRect(30, 340, 700, 100))
        # #
        # self.all_gsa.optimize_group.label_organ = QLabel("Organism name:", self.all_gsa.optimize_group)
        # self.all_gsa.optimize_group.label_organ.setFont(QtGui.QFont('Arial', 12))
        # self.all_gsa.optimize_group.label_organ.setGeometry(QRect(10, 20, 150, 20))
        # self.all_gsa.number_sites_combo = FreeComboBox(self.all_gsa.optimize_group,self.organism_name_list,
        #                                                QRect(215, 20, 120, 20), "Choose Norganism Name",
        #                                                self.change_organism_name)
        # #
        # self.all_gsa.optimize_group.label_codon = QLabel("Codon optimization method:", self.all_gsa.optimize_group)
        # self.all_gsa.optimize_group.label_codon.setFont(QtGui.QFont('Arial', 12))
        # self.all_gsa.optimize_group.label_codon.setGeometry(QRect(10, 50, 200, 20))
        # self.all_gsa.codon_method = FreeComboBox(self.all_gsa.optimize_group,self.codon_methods,
        #                                                QRect(215, 50, 120, 20), "Choose Codon Optimizaion Method",
        #                                                self.change_codon_method)
        # #
        # self.all_gsa.optimize_group.label_gc_content = QLabel("GC Content:", self.all_gsa.optimize_group)
        # self.all_gsa.optimize_group.label_gc_content.setFont(QtGui.QFont('Arial', 12))
        # self.all_gsa.optimize_group.label_gc_content.setGeometry(QRect(370, 20, 150, 20))
        # #
        # self.all_gsa.optimize_group.label_gc_content_min = QLabel("Min:", self.all_gsa.optimize_group)
        # self.all_gsa.optimize_group.label_gc_content_min.setFont(QtGui.QFont('Arial', 12))
        # self.all_gsa.optimize_group.label_gc_content_min.setGeometry(QRect(370, 50, 100, 20))
        # self.all_gsa.gc_min = FreeLineEdit(self.all_gsa.optimize_group,QRect(420, 50, 100, 20), "",str(self.miniGC),self.change_gc_min)
        # #
        # self.all_gsa.optimize_group.label_gc_content_max = QLabel("Max:", self.all_gsa.optimize_group)
        # self.all_gsa.optimize_group.label_gc_content_max.setFont(QtGui.QFont('Arial', 12))
        # self.all_gsa.optimize_group.label_gc_content_max.setGeometry(QRect(540, 50, 100, 20))
        # self.all_gsa.gc_max = FreeLineEdit(self.all_gsa.optimize_group,QRect(590, 50, 100, 20), "",str(self.maxiGC),self.change_gc_max)
        # #
        # self.all_gsa.optimize_group.hide()
        ####


        self.all_gsa.log_label = QLabel("Log:", self.all_gsa)
        self.all_gsa.log_label.setFont(QtGui.QFont('Arial', 12))
        self.all_gsa.log_label.setGeometry(QRect(400, 220, 200, 20))
        self.all_gsa.log = FreeLineEdit(self.all_gsa,QRect(400, 240, 350, 20), "Here you can get some updates","Waiting to start!")
        #
        self.all_gsa.user_guide = PushButton(self.all_gsa,"Open User Guide",QRect(100, 450, 100, 50),"Check Out Our User Guide",self.open_user_guide)
        self.all_gsa.user_guide.button.setFont(QtGui.QFont('Arial', 8))
        self.all_gsa.user_guide.button.setStyleSheet("background-color: turquoise; font-weight: bold")
        #
        self.all_gsa.advanced = PushButton(self.all_gsa,"Optimization Parameters",QRect(300, 295, 150, 30),"Choose advanced options - insert inputh before ",self.optomize_button)
        self.all_gsa.advanced.button.setFont(QtGui.QFont('Arial', 8))
        self.all_gsa.advanced.button.setStyleSheet("background-color: orange; font-weight: bold")
        self.all_gsa.advanced.button.hide()
        #
        self.all_gsa.go = PushButton(self.all_gsa,"GO!",QRect(480, 450, 100, 50)," Start Genomic Stability Analyzer ",self.start)
        self.all_gsa.go.button.setFont(QtGui.QFont('Arial', 12))
        self.all_gsa.go.button.setStyleSheet("background-color: green; font-weight: bold")
        #





if __name__=="__main__":
    print("##########################################################################################\n")
    print("PLEASE WAIT WHILE THE SOFTWARE TURNING ON- IT WILL TAKE A FEW SECONDS\n")
    print("##########################################################################################")
    # Run Program
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec())