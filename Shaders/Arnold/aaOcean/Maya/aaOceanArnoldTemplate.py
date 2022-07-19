import pymel.core as pm
import mtoa.utils as utils
import mtoa.ui.ae.utils as aeUtils
from mtoa.ui.ae.shaderTemplate import ShaderAETemplate

class AEaaOceanTemplate(ShaderAETemplate):
	
	def setup(self):
		self.addSwatch()

		self.beginScrollLayout()

		self.addCustom('message', 'AEshaderTypeNew','AEshaderTypeReplace')

		self.beginLayout("Ocean Parameters", collapse=False)
		self.addControl("resolution", label="Resolution")
		self.addControl("seed", label="Seed")
		self.addControl("oceanScale", label="Ocean Scale")
		self.addControl("time", label="Time (secs)")
		self.addControl("use_uv_input", label="use_uv_input")
		self.addControl("rotateUV", label="rotateUV")
		self.endLayout()
		
		self.beginLayout("Wave Parameters", collapse=False)
		self.addControl("waveHeight", label="Height")
		self.addControl("velocity", label="Size")
		self.addControl("waveSpeed", label="Speed")
		self.addControl("chopAmount", label="Choppiness")
		self.addControl("cutoff", label="Smooth")
		self.endLayout()
		
		self.beginLayout("Wind Parameters", collapse=False)
		self.addControl("windDir", label="Direction")
		self.addControl("damp", label="Reflection")
		self.addControl("windAlign", label="Wind Align")
		self.endLayout()
		
		self.beginLayout("Foam Parameters", collapse=True)
		self.addControl("raw", label="Raw Output")
		self.addControl("invertFoam", label="Invert Foam")
		self.addControl("gamma", label="Gamma")
		self.addControl("brightness", label="Brightness")
		self.addControl("fMin", label="Foam Min")
		self.addControl("fMax", label="Foam Max")
		self.endLayout()
		
		self.beginLayout("File Output", collapse=True)
		self.addControl("writeFile", label="Write Files")
		self.addControl("outputFolder", label="Output Folder")
		self.addControl("postfix", label="Postfix")
		self.addControl("currentFrame", label="Current Frame")
		self.endLayout()


		pm.mel.AEdependNodeTemplate(self.nodeName)