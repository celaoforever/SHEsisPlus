import pygtk
pygtk.require('2.0')
import gtk

class SHEsisGUI:
	def destroy(self, widget, data=None):
		gtk.main_quit()
	
        def AddMenuBar(self):
                self.MenuItems=(
                        ("/File",None,None,0,"<Branch>"),
                        ("/File/Load Case data",None,None,0,None),
                        ("/File/Load Control data",None,None,0,None),
                        ("/File/sep1",None,None,0,"<Separator>"),
                        ("/File/Exit",None,self.destroy,0,None),
                        ("/Help",None,None,0,"<Branch>"),
                        ("/Help/SHEsis Help",None,None,0,None),
                        ("/Help/About SHEsis",None,None,0,None),
                        )
                
                self.AccelGroup=gtk.AccelGroup()
                self.ItemFactory=gtk.ItemFactory(gtk.MenuBar,"<main>",self.AccelGroup)
                self.ItemFactory.create_items(self.MenuItems)
                self.window.add_accel_group(self.AccelGroup)
                self.MenuBar=self.ItemFactory.get_widget("<main>")
                self.vBox=gtk.VBox(False,1)
                self.vBox.set_border_width(1)
                self.vBox.pack_start(self.MenuBar,False,True,0)
                self.MenuBar.show()
        

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
		self.LabelPloidy.set_text("Ploidy:")
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
             
	def __init__(self):
                #create main window
		self.window=gtk.Window(gtk.WINDOW_TOPLEVEL)
		self.window.set_title("SHEsis")
		self.window.connect("destroy",self.destroy)
		self.window.set_border_width(10)
                self.window.set_title("SHEsis")
                self.window.set_resizable(False)
		#add table to manage layout
		self.table=gtk.Table(3,5,False)
		self.table.set_row_spacings(5)
		self.table.set_col_spacings(10)
                #add menu bar
                self.AddMenuBar()
                #add components
                self.AddAnalysisType()
                self.AddPloidy()
                self.AddMarkerName()
		self.AddPhenotype()
		self.AddMask()
                #add table to window and show
		self.vBox.pack_start(self.table,False,True,20)
                self.window.add(self.vBox)
		self.table.show()
                self.vBox.show()
		self.window.show()
	
	def main(self):
		gtk.main()

if __name__ == "__main__":
	gui=SHEsisGUI()
	gui.main()
