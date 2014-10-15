import pygtk
pygtk.require('2.0')
import gtk

class SHEsisGUI:
	def destroy(self, widget, data=None):
		gtk.main_quit()
	
	def __init__(self):
		self.window=gtk.Window(gtk.WINDOW_TOPLEVEL)
		self.window.set_title("SHEsis")
		self.window.connect("destroy",self.destroy)
		self.window.set_border_width(10)
		#add table to manage layout
		self.table=gtk.Table(3,5,False)
		self.table.set_row_spacings(5)
		self.table.set_col_spacings(10)
		#Analysis type
		self.LabelChooseAnalysis=gtk.Label()
		self.LabelChooseAnalysis.set_text("Choose analysis:")
		self.table.attach(self.LabelChooseAnalysis,0,1,0,1)
		self.LabelChooseAnalysis.show()
		self.CheckAssoc=gtk.CheckButton("Association Test")
		self.table.attach(self.CheckAssoc,1,2,0,1)
		self.CheckAssoc.show()
		self.CheckHWE=gtk.CheckButton("Hardy-Weinberg Equilibrium")
		self.table.attach(self.CheckHWE,2,3,0,1)
		self.CheckHWE.show()
		self.CheckLD=gtk.CheckButton("Linkage Disequilibrium")
		self.table.attach(self.CheckLD,3,4,0,1)
		self.CheckLD.show()
		self.CheckHap=gtk.CheckButton("Haplotype")
		self.table.attach(self.CheckHap,4,5,0,1)
		self.CheckHap.show()		
		#ploidy
		self.LabelPloidy=gtk.Label()
		self.LabelPloidy.set_text("Ploidy:")
		self.table.attach(self.LabelPloidy,0,1,1,2)
		self.LabelPloidy.show()
		
		

		self.window.add(self.table)
		self.table.show()
		self.window.show()
	
	def main(self):
		gtk.main()

if __name__ == "__main__":
	gui=SHEsisGUI()
	gui.main()
