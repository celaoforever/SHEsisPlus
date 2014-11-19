import pygtk
pygtk.require('2.0')
import gtk
import os
import platform
import subprocess
import shlex
from subprocess import Popen, PIPE
import glib
import gobject
import threading

class SHEsisGUI:
	def destroy(self, widget, data=None):
		gtk.main_quit()

	def alert(self,msg):
		message=gtk.MessageDialog(None,0,gtk.MESSAGE_ERROR,gtk.BUTTONS_OK,"Error")
		message.set_markup(msg)
		message.run()
		message.destroy()
        
        def read_output(self,view,buffer,cmd):
                args = shlex.split(cmd)
                self.proc = Popen(args,stdout=PIPE,stderr=subprocess.STDOUT)   
                while 1:
                        line=self.proc.stdout.read(1)
                        if not line:
                                break
                        gtk.gdk.threads_enter()
                        iter=buffer.get_end_iter()
                        buffer.insert(iter,line)
                        view.scroll_to_mark(buffer.get_insert(),0.1)
                        gtk.gdk.threads_leave()

        def WriteToFile(self,textbuffer,filename):
                file=open(filename,'w')
                text=textbuffer.get_text(*textbuffer.get_bounds())
                file.write(text)
                file.close()

	def LoadCaseFile(self,widget,data=None):
		chooser = gtk.FileChooserDialog(title="Load case data...",action=gtk.FILE_CHOOSER_ACTION_OPEN,buttons=(gtk.STOCK_CANCEL,gtk.RESPONSE_CANCEL,gtk.STOCK_OPEN,gtk.RESPONSE_OK))
		response = chooser.run()
		if response == gtk.RESPONSE_OK:
			file=open(chooser.get_filename(),'r')
			self.textbufferCase.set_text(file.read())
                        self.casepath=chooser.get_filename()
		chooser.destroy()

	def LoadCtrlFile(self,widget,data=None):
		chooser = gtk.FileChooserDialog(title="Load control data...",action=gtk.FILE_CHOOSER_ACTION_OPEN,buttons=(gtk.STOCK_CANCEL,gtk.RESPONSE_CANCEL,gtk.STOCK_OPEN,gtk.RESPONSE_OK))
		response = chooser.run()
		if response == gtk.RESPONSE_OK:
			file=open(chooser.get_filename(),'r')
			self.textbufferCtrl.set_text(file.read())
			self.ctrlpath=chooser.get_filename()
		chooser.destroy()

	def LoadQTLFile(self,widget,data=None):
		chooser = gtk.FileChooserDialog(title="Load QTL data...",action=gtk.FILE_CHOOSER_ACTION_OPEN,buttons=(gtk.STOCK_CANCEL,gtk.RESPONSE_CANCEL,gtk.STOCK_OPEN,gtk.RESPONSE_OK))
		response = chooser.run()
		if response == gtk.RESPONSE_OK:
			file=open(chooser.get_filename(),'r')
			self.textbufferQTL.set_text(file.read())
			self.qtlpath=chooser.get_filename()
		chooser.destroy()

        
        def AddAnalysisType(self):
  		self.LabelChooseAnalysis=gtk.Label()
		self.LabelChooseAnalysis.set_text("Choose analysis:")
                self.LabelChooseAnalysis.set_alignment(1,0.5)
		self.table.attach(self.LabelChooseAnalysis,0,1,0,1)
		self.LabelChooseAnalysis.show()
		self.CheckAssoc=gtk.CheckButton("Association Test")
                self.CheckAssoc.set_alignment(1,0.5)
		self.table.attach(self.CheckAssoc,1,2,0,1)
		self.CheckAssoc.show()
		self.CheckHWE=gtk.CheckButton("Hardy-Weinberg Equilibrium")
                self.CheckHWE.set_alignment(1,0.5)
		self.table.attach(self.CheckHWE,2,3,0,1)
		self.CheckHWE.show()
		self.CheckLD=gtk.CheckButton("Linkage Disequilibrium")
                self.CheckLD.set_alignment(1,0.5)
		self.table.attach(self.CheckLD,3,4,0,1)
		self.CheckLD.show()
		self.CheckHap=gtk.CheckButton("Haplotype")
                self.CheckHap.set_alignment(1,0.5)
		self.table.attach(self.CheckHap,4,5,0,1)
		self.CheckHap.show()    
          
        def AddPloidy(self):
 		self.LabelPloidy=gtk.Label()
		self.LabelPloidy.set_text("Enter ploidy:")
                self.LabelPloidy.set_alignment(1,0.5)
		self.table.attach(self.LabelPloidy,0,1,1,2)
		self.LabelPloidy.show() 
                adj=gtk.Adjustment(2,1,2048,1,10,0)
                self.SpinPloidy=gtk.SpinButton(adj,0.1,0)
                self.table.attach(self.SpinPloidy,1,2,1,2)
                self.SpinPloidy.show()

        def AddMarkerName(self):
 		self.LabelMarkerName=gtk.Label()
		self.LabelMarkerName.set_text("Marker Name (Optional):")
                self.LabelMarkerName.set_alignment(1,0.5)
		self.table.attach(self.LabelMarkerName,2,3,1,2)
		self.LabelMarkerName.show()                 
                self.TextMarkerName=gtk.Entry()
                self.table.attach(self.TextMarkerName,3,5,1,2)
                self.TextMarkerName.show()

	def ShowCaseCtrlData(self):
		self.LabelDataCase.show()
		self.TextviewDataCase.show()
		self.scrolledwindowCase.show()
		self.LabelDataCtrl.show()
		self.TextviewDataCtrl.show()
		self.scrolledwindowCtrl.show()

	def HideCaseCtrlData(self):
		self.LabelDataCase.hide()
		self.TextviewDataCase.hide()
		self.scrolledwindowCase.hide()
		self.LabelDataCtrl.hide()
		self.TextviewDataCtrl.hide()
		self.scrolledwindowCtrl.hide()
		
	def ShowQTLData(self):
		self.LabelQTL.show()
		self.TextviewQTL.show()
		self.scrolledwindowQTL.show()

	def HideQTLData(self):
		self.LabelQTL.hide()
		self.TextviewQTL.hide()
		self.scrolledwindowQTL.hide()

	def PhenotypeOnChange(self,widget):
		model=self.ComboPhenotype.get_model()
		index=self.ComboPhenotype.get_active()
		if model[index][0] == "Case/Control" :
			self.HideQTLData()
			self.ShowCaseCtrlData()
			self.ComboLD.set_button_sensitivity(gtk.SENSITIVITY_ON)
		elif model[index][0] == "Quantitative Trait":
			self.HideCaseCtrlData()
			self.ShowQTLData()
			self.ComboLD.set_active(0)
			self.ComboLD.set_button_sensitivity(gtk.SENSITIVITY_OFF)
		else:
			raise Exception("Unknown phenotype")	
		
        def AddPhenotype(self):
 		self.LabelPhenotype=gtk.Label()
		self.LabelPhenotype.set_text("Choose Phenotype:")
                self.LabelPhenotype.set_alignment(1,0.5)
		self.table.attach(self.LabelPhenotype,0,1,2,3)     
                self.LabelPhenotype.show()
                self.ComboPhenotype=gtk.combo_box_new_text()
                self.ComboPhenotype.append_text("Case/Control")
                self.ComboPhenotype.append_text("Quantitative Trait")
                self.ComboPhenotype.set_active(0)
		self.ComboPhenotype.connect("changed",self.PhenotypeOnChange)
                self.table.attach(self.ComboPhenotype,1,2,2,3)
                self.ComboPhenotype.show()

        def AddMask(self):
 		self.LabelMask=gtk.Label()
		self.LabelMask.set_text("Mask for hap-analysis (Optional):")
                self.LabelMask.set_alignment(1,0.5)
		self.table.attach(self.LabelMask,2,3,2,3)
                self.LabelMask.show()          
                self.TextMaskName=gtk.Entry()
                self.table.attach(self.TextMaskName,3,5,2,3)
                self.TextMaskName.show() 

	def AddLD(self):
		self.LabelLD=gtk.Label()
		self.LabelLD.set_text("Calculate LD in:")
		self.LabelLD.set_alignment(1,0.5)
		self.table.attach(self.LabelLD,0,1,3,4)
		self.LabelLD.show()
		self.ComboLD=gtk.combo_box_new_text()
		self.ComboLD.append_text("All samples")
		self.ComboLD.append_text("Just cases")
		self.ComboLD.append_text("Just controls")
		self.ComboLD.set_active(0)
		self.table.attach(self.ComboLD,1,2,3,4)
		self.ComboLD.show()

	def AddLFT(self):
		self.LabelLFT=gtk.Label()
		self.LabelLFT.set_text("Lowest freq for hap-analysis:")
		self.LabelLFT.set_alignment(1,0.5)
		self.table.attach(self.LabelLFT,2,3,3,4)
		self.LabelLFT.show()
		self.TextLFT=gtk.Entry()
		self.TextLFT.set_text("0.03")
		self.table.attach(self.TextLFT,3,5,3,4)
		self.TextLFT.show()
	
	def AddOutput(self):
		self.LabelOutput=gtk.Label()
		self.LabelOutput.set_text("Output format:")
		self.LabelOutput.set_alignment(1,0.5)
		self.table.attach(self.LabelOutput,2,3,5,6)
		self.LabelOutput.show()
		self.ComboOutput=gtk.combo_box_new_text()
		self.ComboOutput.append_text("html")
		self.ComboOutput.append_text("txt")
		self.ComboOutput.set_active(0)
		self.table.attach(self.ComboOutput,3,5,5,6)
		self.ComboOutput.show()

	def AddAlgorithm(self):
		self.LabelAlg=gtk.Label()
		self.LabelAlg.set_text("Algorithm for hap-analysis:")
		self.LabelAlg.set_alignment(1,0.5)
		self.table.attach(self.LabelAlg,2,3,4,5)
		self.LabelAlg.show()
		self.ComboAlg=gtk.combo_box_new_text()
		self.ComboAlg.append_text("Expectation maximization")
		self.ComboAlg.append_text("SAT-based")
		self.ComboAlg.set_active(0)
		self.table.attach(self.ComboAlg,3,5,4,5)
		self.ComboAlg.show()
		
	def AddMultiComp(self):
		self.LabelMultiComp=gtk.Label()
		self.LabelMultiComp.set_text("Multiple comparision:")
		self.LabelMultiComp.set_alignment(1,0.5)
		self.table.attach(self.LabelMultiComp,0,1,4,5)
		self.LabelMultiComp.show()
		self.CheckMultiComp=gtk.CheckButton("Sidak,FDR,Holm's correction")
		self.table.attach(self.CheckMultiComp,1,2,4,5)
		self.CheckMultiComp.show()
		
	def AddPermutation(self):
		self.LabelPermutation=gtk.Label()
		self.LabelPermutation.set_text("Number of permutations:")
		self.LabelPermutation.set_alignment(1,0.5)
		self.table.attach(self.LabelPermutation,0,1,5,6)
		self.LabelPermutation.show()
                adj=gtk.Adjustment(0,0,99999999,100,1000,0)
                self.SpinPermutation=gtk.SpinButton(adj,0.1,0)
                self.table.attach(self.SpinPermutation,1,2,5,6)
                self.SpinPermutation.show()


		
	def AddDataInput(self):
		self.TableData=gtk.Table(3,3,False)
		self.TableData.set_col_spacing(1,20)
		self.TableData.set_row_spacing(1,400)
		#case/contrl
		self.LabelDataCase=gtk.Label()
		self.LabelDataCase.set_text("Input case data:")
		self.LabelDataCase.set_alignment(0,0.5)
		self.TableData.attach(self.LabelDataCase,0,1,0,1)
		self.LabelDataCtrl=gtk.Label()
		self.LabelDataCtrl.set_text("Input control data:")
		self.LabelDataCtrl.set_alignment(0,0.5)
		self.TableData.attach(self.LabelDataCtrl,2,3,0,1)
		self.scrolledwindowCase=gtk.ScrolledWindow()
		self.scrolledwindowCase.set_policy(gtk.POLICY_AUTOMATIC,gtk.POLICY_AUTOMATIC)
		self.TextviewDataCase=gtk.TextView()
		self.TextviewDataCase.set_border_window_size(gtk.TEXT_WINDOW_TOP,10)
		self.textbufferCase=self.TextviewDataCase.get_buffer()
		self.scrolledwindowCase.add(self.TextviewDataCase)
		self.TableData.attach(self.scrolledwindowCase,0,1,1,3)
		self.scrolledwindowCtrl=gtk.ScrolledWindow()
		self.scrolledwindowCtrl.set_policy(gtk.POLICY_AUTOMATIC,gtk.POLICY_AUTOMATIC)
		self.TextviewDataCtrl=gtk.TextView()
		self.TextviewDataCtrl.set_border_window_size(gtk.TEXT_WINDOW_TOP,10)
		self.textbufferCtrl=self.TextviewDataCtrl.get_buffer()
		self.scrolledwindowCtrl.add(self.TextviewDataCtrl)
		self.TableData.attach(self.scrolledwindowCtrl,2,3,1,3)		
		#qtl
		self.LabelQTL=gtk.Label()
		self.LabelQTL.set_text("Input QTL data:")
		self.LabelQTL.set_alignment(0,0.5)
		self.TableData.attach(self.LabelQTL,0,1,0,1)
		self.scrolledwindowQTL=gtk.ScrolledWindow()
		self.scrolledwindowQTL.set_policy(gtk.POLICY_AUTOMATIC,gtk.POLICY_AUTOMATIC)
		self.TextviewQTL=gtk.TextView()
		self.TextviewQTL.set_border_window_size(gtk.TEXT_WINDOW_TOP,10)
		self.textbufferQTL=self.TextviewQTL.get_buffer()
		self.scrolledwindowQTL.add(self.TextviewQTL)
		self.TableData.attach(self.scrolledwindowQTL,0,3,1,3)
		#show
		self.table.attach(self.TableData,0,5,6,7)
		self.ShowCaseCtrlData()
		self.TableData.show()

	def ParseArgs(self):
                os=platform.system()
                if os =="Windows":
                        self.argument="SHEsisPlus "
                elif os == "Darwin" or os == "Linux":
                        self.argument="./SHEsis "
                else:
                        raise Exception("Unknown operating system")
		analysisCount=0;
		if self.CheckAssoc.get_active():
			self.argument+=" --assoc"
			analysisCount+=1
		if self.CheckHWE.get_active():
			self.argument+=" --hwe"
			analysisCount+=1
		marker=self.TextMarkerName.get_text()
		ploidy=str(self.SpinPloidy.get_value_as_int())
		self.argument+=" --ploidy "+ploidy
		PhenotypeModel=self.ComboPhenotype.get_model()
		PhenotypeIndex=self.ComboPhenotype.get_active()
		if PhenotypeModel[PhenotypeIndex][0] == "Quantitative Trait":
			self.argument+=" --qtl"
                        if self.qtlpath =="":
                                self.WriteToFile(self.textbufferQTL,"qtl.txt")
                                self.qtlpath="qtl.txt"
                        self.argument+=" --input "+self.qtlpath
                elif PhenotypeModel[PhenotypeIndex][0] == "Case/Control":
                        if self.casepath=="":
                                self.WriteToFile(self.textbufferCase,"case.txt")
                                self.casepath="case.txt"
                        if self.ctrlpath=="":
                                self.WriteToFile(self.textbufferCtrl,"ctrl.txt")
                                self.ctrlpath="ctrl.txt"
                        self.argument+=" --input-case "+self.casepath
                        self.argument+=" --input-ctrl "+self.ctrlpath
                else:
                        raise Exception("Unknown phenotype")
                
		if self.CheckMultiComp.get_active():
			self.argument+=" --adjust"
		permutation=self.SpinPermutation.get_value_as_int()
		if permutation > 0:
			self.argument+=" --permutation "+str(permutation)
		if marker!="":
			self.argument+=' --snpname-line "'+marker+'"'
		if self.CheckLD.get_active():
			analysisCount+=1
			LDModel=self.ComboLD.get_model()
			LDIndex=self.ComboLD.get_active()
			ld=LDModel[LDIndex][0]
			if ld == "All samples":
				self.argument+=" --ld"
			elif ld == "Just cases":
				self.argument+=" --ld-in-case"
			elif ld == "Just controls":
				self.argument+=" --ld-in-ctrl"
			else:
				raise Exception("Unknown ld type")
		if self.CheckHap.get_active():
			analysisCount+=1
			HapModel=self.ComboAlg.get_model()
			HapIndex=self.ComboAlg.get_active()
			hap=HapModel[HapIndex][0]
			if hap == "Expectation maximization":
				self.argument+=" --haplo-EM"
			elif hap == "SAT-based":
				self.argument+=" --haplo-SAT"
			else:
				raise Exception("Unknwon haplotype algorithm")
			mask=self.TextMaskName.get_text()
			if mask!="":
				self.argument+=' --mask "'+mask+'"'
			lft=self.TextLFT.get_text()
			if lft!="":
				self.argument+=" --lft "+lft	
		OutputModel=self.ComboOutput.get_model()
		OutputIndex=self.ComboOutput.get_active()
		Output=OutputModel[OutputIndex][0]
		if Output == "txt":
			self.argument+=" --report-txt"
		if analysisCount ==0:
			self.alert("At leaset one analysis should be selected")
			return -1
		return 0

        def TerminateProcess(self,widget):
                if self.proc.poll() is None:
                        self.proc.terminate()
                        iter=self.textbufferRuntime.get_end_iter()
                        self.textbufferRuntime.insert(iter,"\nTerminated by user")
                        self.TextviewRuntime.scroll_to_mark(self.textbufferRuntime.get_insert(),0.1)
                else:
                        self.SHEsisDialog.destroy()
                
	def CreateDialog(self):
		self.SHEsisDialog=gtk.Dialog("SHEsisPlus",self.window,gtk.DIALOG_DESTROY_WITH_PARENT,None)
		self.LabelRuntime=gtk.Label(" SHEsisPlus runtime output:")
		self.LabelRuntime.set_alignment(0,0.5)
		self.LabelRuntime.show()
		self.SHEsisDialog.vbox.pack_start(self.LabelRuntime,True,True,10)
		self.TableRuntime=gtk.Table(4,3)
		self.TableRuntime.set_row_spacing(0,300)
		self.TableRuntime.set_row_spacing(1,6)
		self.TableRuntime.set_col_spacing(0,200)
		self.TableRuntime.set_col_spacing(1,200)
		self.SHEsisDialog.action_area.pack_start(self.TableRuntime,True,True,5)			
		self.scrolledwindowRuntime=gtk.ScrolledWindow()
		self.scrolledwindowRuntime.set_policy(gtk.POLICY_AUTOMATIC,gtk.POLICY_AUTOMATIC)
		self.TextviewRuntime=gtk.TextView()
		self.TextviewRuntime.set_editable(False)
		self.TextviewRuntime.set_border_window_size(gtk.TEXT_WINDOW_TOP,10)
		self.textbufferRuntime=self.TextviewRuntime.get_buffer()
		self.scrolledwindowRuntime.add(self.TextviewRuntime)
		self.TableRuntime.attach(self.scrolledwindowRuntime,0,3,0,2)
		self.ButtonRuntime=gtk.Button("Terminate")
                self.ButtonRuntime.connect("clicked",self.TerminateProcess)
		self.TableRuntime.attach(self.ButtonRuntime,1,2,3,4)
		self.ButtonRuntime.show()
		self.scrolledwindowRuntime.show()
		self.TextviewRuntime.show()
		self.TableRuntime.show()
		self.SHEsisDialog.show()
	
	def About(self,widget,data=None):
		self.SHEsisAbout=gtk.Dialog("About SHEsisPlus",self.window,gtk.DIALOG_DESTROY_WITH_PARENT,(gtk.STOCK_OK,gtk.RESPONSE_OK))
		self.LabelAbout=gtk.Label("\n\n\t\t\t\t\t\tSHEsisPlus, by Jiawei Shen\n\n\t\tIf you have any problem, please email to jiawei.shen@outlook.com\t\t \n\n\n")
		self.SHEsisAbout.set_resizable(False)
		self.LabelAbout.set_alignment(0,0.5)
		self.LabelAbout.show()
		self.SHEsisAbout.vbox.pack_start(self.LabelAbout,True,True,10)
		r=self.SHEsisAbout.run()
		if r == gtk.RESPONSE_OK:
			self.SHEsisAbout.destroy()
		

	def RunSHEsis(self,widget):
                gobject.threads_init()
                gtk.gdk.threads_init()
		ret=self.ParseArgs()
		if ret == 0:
			self.CreateDialog()
			print self.argument
                        thr=threading.Thread(target=self.read_output,args=(self.TextviewRuntime,self.textbufferRuntime,self.argument))
                        thr.start()

        def reset(self,widget,data=None):
                self.textbufferCase.set_text("")
                self.textbufferCtrl.set_text("")
                self.textbufferQTL.set_text("")
                self.TextMarkerName.set_text("")
                self.TextMaskName.set_text("")
                self.TextLFT.set_text("0.03")
                self.ComboPhenotype.set_active(0)
                self.ComboLD.set_active(0)
                self.ComboAlg.set_active(0)
                self.ComboOutput.set_active(0)
                self.SpinPermutation.set_value(0)
                self.SpinPloidy.set_value(2)
                self.CheckAssoc.set_active(0)
                self.CheckHWE.set_active(0)
                self.CheckLD.set_active(0)
                self.CheckHap.set_active(0)
                self.CheckMultiComp.set_active(0)
	
	def AddCalculate(self):
		self.TableButton=gtk.Table(1,9,False)
		self.TableButton.set_col_spacing(1,310)
		self.TableButton.set_col_spacing(7,290)
		self.TableButton.set_col_spacing(6,40)
		self.ButtonCal=gtk.Button("Calculate")
		self.ButtonCal.set_border_width(5)
		self.ButtonCal.connect("clicked",self.RunSHEsis)
		self.TableButton.attach(self.ButtonCal,3,4,0,1)
		self.ButtonCal.show()
		self.ButtonReset=gtk.Button("  Reset  ")
		self.ButtonReset.set_border_width(5)
                self.ButtonReset.connect("clicked",self.reset)
		self.TableButton.attach(self.ButtonReset,4,5,0,1)
		self.ButtonReset.show()		
		self.vBox.pack_start(self.TableButton,False,False,0)
		self.TableButton.show()	



        def AddMenuBar(self):
                self.MenuItems=(
                        ("/File",None,None,0,"<Branch>"),
                        ("/File/Load Case data",None,self.LoadCaseFile,0,None),
                        ("/File/Load Control data",None,self.LoadCtrlFile,0,None),
			("/File/Load QTL data",None,self.LoadQTLFile,0,None),
                        ("/File/sep1",None,None,0,"<Separator>"),
                        ("/File/Exit",None,self.destroy,0,None),
			("/Edit",None,None,0,"<Branch>"),
			("/Edit/Reset",None,self.reset,0,None),
                        ("/Help",None,None,0,"<Branch>"),
                        ("/Help/SHEsisPlus Help",None,None,0,None),
                        ("/Help/About SHEsisPlus",None,self.About,0,None),
                        )
                
                self.AccelGroup=gtk.AccelGroup()
                self.ItemFactory=gtk.ItemFactory(gtk.MenuBar,"<main>",self.AccelGroup)
                self.ItemFactory.create_items(self.MenuItems)
                self.window.add_accel_group(self.AccelGroup)
                self.MenuBar=self.ItemFactory.get_widget("<main>")
		self.vBox.pack_start(self.MenuBar,False,True,0)
                self.MenuBar.show()

		
	def __init__(self):
		self.casepath=""
		self.ctrlpath=""
		self.qtlpath=""
                #create main window
		self.window=gtk.Window(gtk.WINDOW_TOPLEVEL)
		self.window.set_title("SHEsisPlus")
		self.window.connect("destroy",self.destroy)
		self.window.set_border_width(10)
                self.window.set_title("SHEsisPlus")
                self.window.set_resizable(False)
		#create vbox
		self.vBox=gtk.VBox(False,1)
                self.vBox.set_border_width(1)
                #add menu bar
                self.AddMenuBar()
		#add table to manage layout
		self.table=gtk.Table(3,5,False)
		self.table.set_row_spacings(5)
		self.table.set_col_spacings(10)
		self.vBox.pack_start(self.table,False,True,10)
                #add components
                self.AddAnalysisType()
                self.AddPloidy()
                self.AddMarkerName()
		self.AddPhenotype()
		self.AddMask()
		self.AddLD()
		self.AddLFT()
		self.AddOutput()
		self.AddAlgorithm()
		self.AddMultiComp()
		self.AddPermutation()
		self.AddDataInput()
		self.AddCalculate()
                #add table to window and show		
                self.window.add(self.vBox)
		self.table.show()
                self.vBox.show()
		self.window.show()

	def main(self):
		gtk.main()

if __name__ == "__main__":
	gui=SHEsisGUI()
	gui.main()
