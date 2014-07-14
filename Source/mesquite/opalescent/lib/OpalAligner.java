/* Opalescent source code.  Copyright 2012 Travis Wheeler and David Maddison.
Version 2.10, January 2012.

Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
Mesquite's web site is http://mesquiteproject.org

This source code and its compiled class files are free and modifiable under the terms of 
GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)
*/
package mesquite.opalescent.lib;


import opal.makers.AlignmentMaker_SingleSequences;
import opal.polish.Polisher;
import opal.tree.Tree;
import opal.tree.Tree.DistanceType;
import opal.align.Aligner;
import opal.align.Alignment;
import opal.align.Aligner.AlignmentType;
import opal.IO.ArgumentHandler;
import opal.exceptions.GenericOpalException;
import opal.IO.CostMatrix;
import opal.IO.OpalLogWriter;
import opal.IO.SequenceConverter;
import com.traviswheeler.libs.Logger;
import mesquite.lib.*;


public class OpalAligner {
	public static boolean defaultFastAlgorithm = true;
	public static boolean defaultPolish = true;
	public static boolean defaultVerbose = true;
	boolean fastAlgorithm = defaultFastAlgorithm;
	boolean polish = defaultPolish;
	boolean isProtein = true;
	boolean verbose = defaultVerbose;
	Logger logger=null;
	String costName = DEFAULTPROTCOSTNAME;
	public static String DEFAULTPROTCOSTNAME = "BLOSUM62";
	String protCostName = "BLOSUM62";
	String extraArguments = "";
	
	int gamma, gammaTerm, lambda, lambdaTerm;

	int protGamma, protGammaTerm, protLambda, protLambdaTerm;
	int dnaGamma, dnaGammaTerm, dnaLambda, dnaLambdaTerm;

	int dnaSubAG, dnaSubCT;
	static int dnaTransversion=100;

	public OpalAligner() {
		setDefaultCosts();
		initCosts();
	}
	
	public OpalAligner(boolean fastAlgorithm, boolean polish, Logger logger) {
		this.fastAlgorithm = fastAlgorithm;
		this.polish = polish;
		this.logger = logger;
		setDefaultCosts();
		initCosts();
	}
	

	public String getOpalVersion() {
		OpalLogWriter OLW = new OpalLogWriter();
		if (OLW!=null) {
			String s = OLW.getClass().getPackage().getImplementationVersion();
			OLW=null;
			return s;
		}
		return "";
	}
	
	public int[][] getAlignment (int[][] input) {
		OpalLogWriter.setLogger(logger);
		CostMatrix.isDNA = !isProtein;

		initCosts();
		
		int[][] ret = null;
		
		try {
			logger.logln("Opal: initialization.");
			CostMatrix.initialize(costName, gamma, gammaTerm, lambda, lambdaTerm);
	
			if (fastAlgorithm) {
				Aligner.alignmentMethod = AlignmentType.profile;
			} else {
				Aligner.alignmentMethod = AlignmentType.mixed;
			}
			
			if (!polish) {
				Polisher.polishMethod = null;
			}
			CostMatrix.dnaSubAG = dnaSubAG;
			CostMatrix.dnaSubCT = dnaSubCT;

			
			Alignment.setAlphabetLength(CostMatrix.getChars().length);
			SequenceConverter sc = new SequenceConverter(CostMatrix.getChars());
			Aligner.setSequenceConverter(sc);		
	
			Aligner.setParams();
			
	// now add in any extraArguments.  This should be done AFTER all other options are fed into Opal
			if (StringUtil.notEmpty(extraArguments)) {
				ArgumentHandler argHandler = new ArgumentHandler(extraArguments);
				// Not doing anything with this object. 
				// There are some getter functions that return a subset of all the
				// values that were passed in the string, but these happen to be ones
				// that I think should be overridden by setting directly assigned in Mesquite anyway.
				/// ... we can discuss
				
			}

			// only set these if they weren't set with custom arguments
			if (Tree.distanceType == null)
				Tree.distanceType = DistanceType.kmer_normcost;    
			if (Tree.iterations == -1)
				Tree.iterations = 2; 

			
	//if aligning alignments use AlignmentMaker_TwoAlignments, and pass two arguments to am.initialize
			AlignmentMaker_SingleSequences am = new AlignmentMaker_SingleSequences();
			AlignmentMaker_SingleSequences.writeResult = false;
			
			int[][] seqs = Alignment.getDegappedSeqs(input);
			
			am.initialize(seqs, null);
			logger.logln("Opal: building alignment.");
			//ret = am.buildAlignment(costName, !verbose /*quiet*/, false /*toUpper*/);
			ret = am.buildAlignment(costName, 1 /*standard verbosity*/, false /*toUpper*/);
			
			
			
		} catch (GenericOpalException e) {
			logger.logln("Opal failed to complete alignment. See log for details");
		} catch (Exception e) {
			logger.logln("Error running Opal. (" + e.getMessage() + ")");
		}
		
		return ret; 
	}
	
	public boolean getFastAlgorithm() {
		return fastAlgorithm;
	}
	public void setFastAlgorithm(boolean fastAlgorithm) {
		this.fastAlgorithm = fastAlgorithm;
	}
	public boolean getPolish() {
		return polish;
	}
	public void setPolish(boolean polish) {
		this.polish = polish;
	}
	public boolean isProtein() {
		return isProtein;
	}
	public void setProtein(boolean isProtein) {
		this.isProtein = isProtein;
		initCosts();
	}

	public Logger getLogger() {
		return logger;
	}

	public void setLogger(Logger logger) {
		this.logger = logger;
	}

	public String getCostName() {
		return costName;
	}

	public void setCostName(String costName) {
		this.costName = costName;
	}

	public boolean isVerbose() {
		return verbose;
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

	public void initCosts(){
		if (isProtein) {
			gamma = protGamma;
			gammaTerm = protGammaTerm;
			lambda = protLambda;
			lambdaTerm = protLambdaTerm;   
			costName = protCostName;
		} else {
			gamma = dnaGamma;
			gammaTerm = dnaGammaTerm;
			lambda = dnaLambda;   
			lambdaTerm = dnaLambdaTerm;   
			costName = "DNA";
		}		
	}
	public void setToDefaults() {
		setDefaultCosts();
		verbose = defaultVerbose;
		fastAlgorithm = defaultFastAlgorithm;
		polish = defaultPolish;

	}
	public void setDefaultCosts(){
			protGamma = CostMatrix.protDefaultGamma;
			protGammaTerm = CostMatrix.protDefaultGammaTerm;
			protLambda = CostMatrix.protDefaultLambda;
			protLambdaTerm = CostMatrix.protDefaultLambdaTerm;   
			dnaGamma = CostMatrix.dnaDefaultGamma;
			dnaGammaTerm = CostMatrix.dnaDefaultGammaTerm;
			dnaLambda = CostMatrix.dnaDefaultLambda;   
			dnaLambdaTerm = CostMatrix.dnaDefaultLambdaTerm;   
			protCostName = DEFAULTPROTCOSTNAME;
			dnaSubAG = CostMatrix.dnaSubAG;
			dnaSubCT = CostMatrix.dnaSubCT;
	}
	
	public void setProtCosts(int protGamma, int protGammaTerm, int protLambda, int protLambdaTerm){
		this.protGamma = protGamma;
		this.protGammaTerm = protGammaTerm;
		this.protLambda = protLambda;
		this.protLambdaTerm = protLambdaTerm;
	}
	public void setdnaCosts(int dnaGamma, int dnaGammaTerm, int dnaLambda, int dnaLambdaTerm){
		this.dnaGamma = dnaGamma;
		this.dnaGammaTerm = dnaGammaTerm;
		this.dnaLambda = dnaLambda;
		this.dnaLambdaTerm = dnaLambdaTerm;
	}

	public int getGamma() {
		return gamma;
	}

	public int getGammaTerm() {
		return gammaTerm;
	}

	public int getLambda() {
		return lambda;
	}

	public int getLambdaTerm() {
		return lambdaTerm;
	}

	public int getProtGamma() {
		return protGamma;
	}

	public int getProtGammaTerm() {
		return protGammaTerm;
	}

	public int getProtLambda() {
		return protLambda;
	}

	public int getProtLambdaTerm() {
		return protLambdaTerm;
	}

	public int getDnaGamma() {
		return dnaGamma;
	}

	public int getDnaGammaTerm() {
		return dnaGammaTerm;
	}

	public int getDnaLambda() {
		return dnaLambda;
	}

	public int getDnaLambdaTerm() {
		return dnaLambdaTerm;
	}

	public String getExtraArguments() {
		return extraArguments;
	}

	public void setExtraArguments(String extraArguments) {
		this.extraArguments = extraArguments;
	}

	public String getProtCostName() {
		return protCostName;
	}

	public void setProtCostName(String protCostName) {
		this.protCostName = protCostName;
	}

	public int getDnaSubAG() {
		return dnaSubAG;
	}

	public void setDnaSubAG(int dnaSubAG) {
		this.dnaSubAG = dnaSubAG;
	}

	public int getDnaSubCT() {
		return dnaSubCT;
	}

	public void setDnaSubCT(int dnaSubCT) {
		this.dnaSubCT = dnaSubCT;
	}

	public int getDnaTransversion() {
		return dnaTransversion;
	}


}
