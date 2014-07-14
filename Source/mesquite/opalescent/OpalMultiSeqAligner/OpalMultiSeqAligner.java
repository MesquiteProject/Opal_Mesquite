/* Opalescent source code.  Copyright 2012 Travis Wheeler and David Maddison.
Version 2.10, January 2012.

Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
Mesquite's web site is http://mesquiteproject.org

This source code and its compiled class files are free and modifiable under the terms of 
GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)
 */
package mesquite.opalescent.OpalMultiSeqAligner;

import java.awt.*;
import java.awt.event.*;

import opal.IO.CostMatrix;
import mesquite.opalescent.lib.OpalAligner;
import mesquite.opalescent.lib.OpalUtil;
import mesquite.lib.*;
import mesquite.align.lib.*;


public class OpalMultiSeqAligner extends MultipleSequenceAligner implements com.traviswheeler.libs.Logger, ActionListener {
	boolean fastAlgorithm = true;
	boolean polish = true;
	boolean verbose = true;
	int protGamma, protGammaTerm, protLambda, protLambdaTerm;
	int dnaGamma, dnaGammaTerm, dnaLambda, dnaLambdaTerm;
	int dnaSubAG, dnaSubCT, dnaTransversion;

	String protCostName = OpalAligner.DEFAULTPROTCOSTNAME;

	String extraArguments = "";


	public boolean startJob(String arguments, Object condition, boolean hiredByName) {

		return true;
	}
	public boolean permitSeparateThread(){
		return false;
	}
	private void getDefaultSettings() {
		fastAlgorithm = OpalAligner.defaultFastAlgorithm;
		polish = OpalAligner.defaultPolish;
		verbose = OpalAligner.defaultVerbose;
		protGamma = CostMatrix.protDefaultGamma;
		protGammaTerm = CostMatrix.protDefaultGammaTerm;
		protLambda = CostMatrix.protDefaultLambda;
		protLambdaTerm = CostMatrix.protDefaultLambdaTerm;   
		dnaGamma = CostMatrix.dnaDefaultGamma;
		dnaGammaTerm = CostMatrix.dnaDefaultGammaTerm;
		dnaLambda = CostMatrix.dnaDefaultLambda;   
		dnaLambdaTerm = CostMatrix.dnaDefaultLambdaTerm;   
		protCostName = OpalAligner.DEFAULTPROTCOSTNAME;
		dnaSubAG = CostMatrix.dnaSubAG;
		dnaSubCT = CostMatrix.dnaSubCT;
		dnaTransversion= 100;
	}

	private void getInitialSettings(OpalAligner opalAligner) {
		fastAlgorithm = opalAligner.getFastAlgorithm();
		polish = opalAligner.getPolish();
		verbose = opalAligner.isVerbose();
		protGamma = opalAligner.getProtGamma();
		protGammaTerm = opalAligner.getProtGammaTerm();
		protLambda = opalAligner.getProtLambda();
		protLambdaTerm = opalAligner.getProtLambdaTerm();
		dnaGamma = opalAligner.getDnaGamma();
		dnaGammaTerm = opalAligner.getDnaGammaTerm();
		dnaLambda = opalAligner.getDnaLambda();
		dnaLambdaTerm = opalAligner.getDnaLambdaTerm();
		dnaSubAG = opalAligner.getDnaSubAG();
		dnaSubCT = opalAligner.getDnaSubCT();
		dnaTransversion= opalAligner.getDnaTransversion();
	}
	/*.................................................................................................................*/
	public void processSingleXMLPreference (String tag, String content) {
		if ("fastAlgorithm".equalsIgnoreCase(tag)) {
			fastAlgorithm = MesquiteBoolean.fromTrueFalseString(content);
		} 
		else if ("polish".equalsIgnoreCase(tag)) {
			polish = MesquiteBoolean.fromTrueFalseString(content);
		} 
		else if ("verbose".equalsIgnoreCase(tag)) {
			verbose = MesquiteBoolean.fromTrueFalseString(content);
		} 
		else if ("protGamma".equalsIgnoreCase(tag)) {
			protGamma = MesquiteInteger.fromString(content);
		} 
		else if ("protGammaTerm".equalsIgnoreCase(tag)) {
			protGammaTerm = MesquiteInteger.fromString(content);
		} 
		else if ("protLambda".equalsIgnoreCase(tag)) {
			protLambda = MesquiteInteger.fromString(content);
		} 
		else if ("protLambdaTerm".equalsIgnoreCase(tag)) {
			protLambdaTerm = MesquiteInteger.fromString(content);
		} 
		else if ("dnaGamma".equalsIgnoreCase(tag)) {
			dnaGamma = MesquiteInteger.fromString(content);
		} 
		else if ("dnaGammaTerm".equalsIgnoreCase(tag)) {
			dnaGammaTerm = MesquiteInteger.fromString(content);
		} 
		else if ("dnaLambda".equalsIgnoreCase(tag)) {
			dnaLambda = MesquiteInteger.fromString(content);
		} 
		else if ("dnaLambdaTerm".equalsIgnoreCase(tag)) {
			dnaLambdaTerm = MesquiteInteger.fromString(content);
		} else if ("extraArguments".equalsIgnoreCase(tag)) {
			extraArguments = StringUtil.cleanXMLEscapeCharacters(content);

		} else if ("protCostName".equalsIgnoreCase(tag)) {
			protCostName = StringUtil.cleanXMLEscapeCharacters(content);
		}
		else if ("dnaSubAG".equalsIgnoreCase(tag)) {
			dnaSubAG = MesquiteInteger.fromString(content);
		} 
		else if ("dnaSubCT".equalsIgnoreCase(tag)) {
			dnaSubCT = MesquiteInteger.fromString(content);
		} 


	}

	/*.................................................................................................................*/
	public String preparePreferencesForXML () {
		StringBuffer buffer = new StringBuffer(20);
		StringUtil.appendXMLTag(buffer, 2, "fastAlgorithm", fastAlgorithm);  
		StringUtil.appendXMLTag(buffer, 2, "polish", polish);  
		StringUtil.appendXMLTag(buffer, 2, "verbose", verbose);  

		StringUtil.appendXMLTag(buffer, 2, "protGamma", protGamma);  
		StringUtil.appendXMLTag(buffer, 2, "protGammaTerm", protGammaTerm);  
		StringUtil.appendXMLTag(buffer, 2, "protLambda", protLambda);  
		StringUtil.appendXMLTag(buffer, 2, "protLambdaTerm", protLambdaTerm);  
		StringUtil.appendXMLTag(buffer, 2, "dnaGamma", dnaGamma);  
		StringUtil.appendXMLTag(buffer, 2, "dnaGammaTerm", dnaGammaTerm);  
		StringUtil.appendXMLTag(buffer, 2, "dnaLambda", dnaLambda);  
		StringUtil.appendXMLTag(buffer, 2, "dnaLambdaTerm", dnaLambdaTerm);  

		StringUtil.appendXMLTag(buffer, 2, "extraArguments", extraArguments);  
		StringUtil.appendXMLTag(buffer, 2, "protCostName", protCostName);  

		StringUtil.appendXMLTag(buffer, 2, "dnaSubAG", dnaSubAG);  
		StringUtil.appendXMLTag(buffer, 2, "dnaSubCT", dnaSubCT);  

		return buffer.toString();
	}

	/*.................................................................................................................*/

	protected long[][] alignSequences(long[][] sequences) {
		if (sequences==null)
			return null;
		OpalAligner opalAligner = new OpalAligner(fastAlgorithm, polish, this);  //create the OPAL object
		if (opalAligner==null)
			return null;
		getInitialSettings(opalAligner);
		loadPreferences();


		if (okToInteractWithUser(CAN_PROCEED_ANYWAY, "Querying about options")) //need to check if can proceed
			if (!queryOptions())
				return null;

		logln("Formatting data for Opal...");
		int[] mesquiteRows = new int[sequences[0].length];
		int[][] opalArray = OpalUtil.getOpalArray(sequences, isProtein, mesquiteRows);  //convert to OPAL format




		logln("Creating the Opal Aligner...");
		opalAligner.setVerbose(verbose);
		opalAligner.setdnaCosts(dnaGamma, dnaGammaTerm, dnaLambda, dnaLambdaTerm);
		opalAligner.setProtCosts(protGamma, protGammaTerm, protLambda, protLambdaTerm);
		opalAligner.setProtein(isProtein);
		opalAligner.setExtraArguments(extraArguments);
		opalAligner.setProtCostName(protCostName);
		opalAligner.setPolish(polish);
		opalAligner.setDnaSubAG(dnaSubAG);
		opalAligner.setDnaSubCT(dnaSubCT);


		logln("Asking Opal (version " + opalAligner.getOpalVersion() + ") to align the data...");
		int[][] opalResults = opalAligner.getAlignment(opalArray);  //do the alignment!

		if (opalResults==null)
			return null;
		
		if (opalResults.length != sequences[0].length)
			opalResults = OpalUtil.addBackEmptyRows(opalResults, mesquiteRows );


		logln("Converting the Opal results back to Mesquite's data format...");
		return OpalUtil.getMesquiteArray(opalResults, isProtein);  //convert back to Mesquite format and return it.
	}

	/*.................................................................................................................*/
	/** returns whether the module is compatible with the current OS, Mesquite system, Java VM, and so on.  If false then the module will not be loaded as a possibility at startup*/
	public boolean  compatibleWithSystem(){ 
		try {
			OpalAligner opalAligner = new OpalAligner(false, false, this);  //create the OPAL object
			return true;
		} catch (Exception e) {
			return false;
		}
	}
	Choice protCostChoice;
	Checkbox fastCheckBox;
	Checkbox polishBox;
	SingleLineTextField extraArgumentsField;
	Checkbox verboseBox;
	SimpleIntegerField[][] gapCostFields;
	SimpleIntegerField[][] dnaSubCostsField;

	/*.................................................................................................................*/
	public boolean queryOptions() {
		MesquiteInteger buttonPressed = new MesquiteInteger(1);
		ExtensibleDialog dialog = new ExtensibleDialog(containerOfModule(), "Opal Options",buttonPressed);  //MesquiteTrunk.mesquiteTrunk.containerOfModule()
		if (isProtein) 
			dialog.addLabel("Opal Options for Protein Data");
		else
			dialog.addLabel("Opal Options for Nucleotide Data");
		String helpString = "This module will align sequences using Opal.<br><br>";
		helpString+=" If you publish results based upon an alignment made using this package, please cite the paper that describes Opal: <br><br>";
		helpString+=" Wheeler, T. J., & Kececioglu, J. D. (2007). Multiple alignments by aligning alignments. <i>Bioinformatics</i>, 23, i559-i568.";


		dialog.appendToHelpString(helpString);
		dialog.setHelpURL (this, "../aOpalescentIntro/index.html");

		int[][] dnaSubCosts =new int[1][3];
		dnaSubCosts[0][0] = dnaSubAG;
		dnaSubCosts[0][1] = dnaSubCT;
		dnaSubCosts[0][2] = dnaTransversion;


		int[][] costs =new int[4][1];
		if (!isProtein) {
			costs[0][0] = dnaGamma;
			costs[1][0] = dnaGammaTerm;
			costs[2][0] = dnaLambda;
			costs[3][0] = dnaLambdaTerm;
		} else {
			costs[0][0] = protGamma;
			costs[1][0] = protGammaTerm;
			costs[2][0] = protLambda;
			costs[3][0] = protLambdaTerm;
		}
		dialog.addBlankLine();

		gapCostFields = dialog.addIntegerFieldsMatrix(4, 1, new String[] {"Gap costs"}, new String[] {"Open", "Terminal Open", "Extension", "Terminal Extension"}, new int[] {4, 4, 4, 4}, costs, 1,Integer.MAX_VALUE);
		dialog.addBlankLine();

		if (!isProtein) {
			dnaSubCostsField = dialog.addIntegerFieldsMatrix(1,3, new String[] {"A<->G", "C<->T", "Transversion"}, new String[] {"Sub. Cost"}, new int[] {4}, dnaSubCosts, 1,Integer.MAX_VALUE);
			dnaSubCostsField[0][2].setEnabled(false);
			dialog.addBlankLine();
		} else {
			protCostChoice = dialog.addPopUpMenu("Protein cost matrix", new String[] {"BLOSUM62", "BLOSUM50"}, 0);
			protCostChoice.select(protCostName);
		}

		fastCheckBox = dialog.addCheckBox("fast alignment (recommended for >40 sequences)", fastAlgorithm);
		polishBox = dialog.addCheckBox("polish", polish);
		extraArgumentsField = dialog.addTextField("Other Options:", extraArguments, 44, true);

		verboseBox = dialog.addCheckBox("report progress in more detail", verbose);
		dialog.addBlankLine();
		final Button defaultsButton = dialog.addAListenedButton("Reset to defaults",null, this);
		defaultsButton.setActionCommand("defaults");



		dialog.completeAndShowDialog(true);
		if (buttonPressed.getValue()==0)  {
			fastAlgorithm = fastCheckBox.getState();
			polish = polishBox.getState();
			verbose = verboseBox.getState();
			if (!isProtein) {
				dnaGamma = gapCostFields[0][0].getValue();
				dnaGammaTerm = gapCostFields[1][0].getValue();
				dnaLambda = gapCostFields[2][0].getValue();
				dnaLambdaTerm = gapCostFields[3][0].getValue();
				dnaSubAG = dnaSubCostsField[0][0].getValue();
				dnaSubCT = dnaSubCostsField[0][1].getValue();
			} else {
				protGamma = gapCostFields[0][0].getValue();
				protGammaTerm = gapCostFields[1][0].getValue();
				protLambda = gapCostFields[2][0].getValue();
				protLambdaTerm = gapCostFields[3][0].getValue();
				protCostName = protCostChoice.getSelectedItem();
			}
			extraArguments = extraArgumentsField.getText();



			storePreferences();
		}
		dialog.dispose();
		return (buttonPressed.getValue()==0) ;
	}

	/*.................................................................................................................*/
	/** returns whether this module is requesting to appear as a primary choice */
	public boolean requestPrimaryChoice(){
		return true;  
	}

	public String getName() {
		return "Opal";
	}

	public String getNameForMenuItem() {
		return "Opal...";
	}
	public boolean isPrerelease(){
		return false;
	}


	/*.................................................................................................................*/
	public  void actionPerformed(ActionEvent e) {
		if (e.getActionCommand().equalsIgnoreCase("defaults")) {
			getDefaultSettings();
			if (protCostChoice!=null)
				protCostChoice.select(protCostName);
			if (fastCheckBox!=null)
				fastCheckBox.setState(fastAlgorithm);
			if (polishBox!=null)
				polishBox.setState(polish);
			if (verboseBox!=null) {
				verboseBox.setState(verbose);
			}
			if (extraArgumentsField!=null)
				extraArgumentsField.setText("");
			if (!isProtein) {
				gapCostFields[0][0].setText(""+dnaGamma);
				gapCostFields[1][0].setText(""+dnaGammaTerm);
				gapCostFields[2][0].setText(""+dnaLambda);
				gapCostFields[3][0].setText(""+dnaLambdaTerm);
				dnaSubCostsField[0][0].setText(""+dnaSubAG);
				dnaSubCostsField[0][1].setText(""+dnaSubCT);
				dnaSubCostsField[0][2].setText(""+dnaTransversion);
			} else {
				gapCostFields[0][0].setText(""+protGamma);
				gapCostFields[1][0].setText(""+protGammaTerm);
				gapCostFields[2][0].setText(""+protLambda);
				gapCostFields[3][0].setText(""+protLambdaTerm);
			}


		}
	}

}
