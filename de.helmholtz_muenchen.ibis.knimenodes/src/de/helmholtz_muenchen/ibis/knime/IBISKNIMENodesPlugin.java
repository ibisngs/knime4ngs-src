/**
 *  Copyright (C) 2016 the Knime4NGS contributors.
 *  Website: http://ibisngs.github.io/knime4ngs
 *  
 *  This file is part of the KNIME4NGS KNIME extension.
 *  
 *  The KNIME4NGS extension is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package de.helmholtz_muenchen.ibis.knime;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;

import javax.swing.JOptionPane;

import org.eclipse.jface.preference.IPreferenceStore;
import org.eclipse.swt.widgets.Button;
import org.eclipse.swt.widgets.Display;
import org.eclipse.swt.widgets.Table;
import org.eclipse.swt.widgets.TableItem;
import org.eclipse.ui.plugin.AbstractUIPlugin;
import org.knime.workbench.ui.startup.StartupMessage;
import org.osgi.framework.BundleContext;

import de.helmholtz_muenchen.ibis.utils.BinaryHandler;

/**
 * This is the OSGI bundle activator.
 * 
 * @author hastreiter
 */
public class IBISKNIMENodesPlugin extends AbstractUIPlugin {
	
//	fields of preference page
	public static final String USE_HTE = "hte";
	public static final String DB_FILE = "db_file";
	public static final String THRESHOLD = "threshold";
	public static final String NOTIFY = "notify";
	public static final String EMAIL_HOST = "host";
	public static final String EMAIL_SENDER = "sender";
	public static final String EMAIL_RECEIVER = "receiver";
	public static final String REF_GENOME = "reference sequence";
	public static final String RES_HAPMAP = "HapMap";
	public static final String RES_OMNI = "Omni";
	public static final String RES_1000G_SNPS = "1000G SNPs";
	public static final String RES_1000G_INDELS ="1000G Indels";
	public static final String RES_DBSNP = "dbSNP";
	public static final String RES_MILLS = "Mills";
	
//	default values
	public static final boolean HTE_DEFAULT = false;
	public static final int THRESHOLD_DEFAULT = 1;
	public static final boolean NOTIFY_DEFAULT = false;
//	public static final String EMAIL_HOST_DEFAULT = "outmail.helmholtz-muenchen.de";
//	public static final String EMAIL_SENDER_DEFAULT = "ibis.knime@helmholtz-muenchen.de";
	public static final String EMAIL_HOST_DEFAULT = "";
	public static final String EMAIL_SENDER_DEFAULT = "";
	public static final String EMAIL_RECEIVER_DEFAULT = "";
	
//	tool binaries
	public static final String BCFTOOLS = "bcftools";
	public static final String BOWTIE2 = "bowtie2";
	public static final String BWA = "bwa";
	public static final String FEATURE_COUNTS = "featureCounts";
	public static final String GATK = "GenomeAnalysisTK.jar";
	public static final String KGGSeq = "kggseq.jar";
	public static final String PICARD = "picard.jar";
	public static final String PINDEL = "pindel";
	public static final String PINDEL2VCF = "pindel2vcf";
	public static final String SAMTOOLS = "samtools";
	public static final String SEGEMEHL = "segemehl.x";
	public static final String SNPEFF = "snpEff.jar";
	public static final String SNPSIFT = "SnpSift.jar";
	public static final String STAR = "STAR";
	public static final String VCFTOOLS = "vcftools";
	public static final String VEP = "variant_effect_predictor.pl";
	public static final String VEP_FILTER = "filter_vep.pl";
	
//	tool dependencies, neither in FIELDS nor in PATHS
	public static final String BOWTIE2_BUILD = "bowtie2-build";
	public static final String BOWTIE2_BUILD_S = "bowtie2-build-s";
	public static final String BOWTIE2_BUILD_L = "bowtie2-build-l";
	public static final String BOWTIE2_ALIGN_S = "bowtie2-align-s";
	public static final String BOWTIE2_ALIGN_L = "bowtie2-align-l";
	
//	excluded binaries: integrate binaries by adding them to PATHS and TOOLS in the KNIMEPreferencePage
//	public static final String BFAST = "bfast";
	
	public static String [] FIELDS = {USE_HTE, THRESHOLD, DB_FILE, NOTIFY, EMAIL_HOST, EMAIL_SENDER, EMAIL_RECEIVER};
	public static String [] PATHS = {
			REF_GENOME, RES_HAPMAP, RES_OMNI, RES_1000G_SNPS, RES_1000G_INDELS , RES_DBSNP, RES_MILLS, 
			DB_FILE,
			BCFTOOLS, BOWTIE2, BWA,
			FEATURE_COUNTS,
			GATK,
			KGGSeq,
			PICARD, PINDEL, PINDEL2VCF,
			SAMTOOLS, SEGEMEHL, SNPEFF, SNPSIFT, STAR,
			VCFTOOLS, VEP, VEP_FILTER};
	
	public static final HashMap<String, Boolean> TOOLS;
	static {
		TOOLS = new HashMap<>();
		TOOLS.put(IBISKNIMENodesPlugin.BCFTOOLS, true);
		TOOLS.put(IBISKNIMENodesPlugin.BOWTIE2,true);
		TOOLS.put(IBISKNIMENodesPlugin.BWA,true);
		TOOLS.put(IBISKNIMENodesPlugin.FEATURE_COUNTS,true);
		TOOLS.put(IBISKNIMENodesPlugin.GATK,false);
		TOOLS.put(IBISKNIMENodesPlugin.KGGSeq,false);
		TOOLS.put(IBISKNIMENodesPlugin.PICARD,false);
		TOOLS.put(IBISKNIMENodesPlugin.PINDEL,true);
		TOOLS.put(IBISKNIMENodesPlugin.PINDEL2VCF,true);
		TOOLS.put(IBISKNIMENodesPlugin.SAMTOOLS,true);
		TOOLS.put(IBISKNIMENodesPlugin.SEGEMEHL,true);
		TOOLS.put(IBISKNIMENodesPlugin.SNPEFF, false);
		TOOLS.put(IBISKNIMENodesPlugin.SNPSIFT, false);
		TOOLS.put(IBISKNIMENodesPlugin.STAR,true);
		TOOLS.put(IBISKNIMENodesPlugin.VCFTOOLS,false);
		TOOLS.put(IBISKNIMENodesPlugin.VEP,false);
		TOOLS.put(IBISKNIMENodesPlugin.VEP_FILTER, false);
	}
	
	/**
     * The shared instance.
     */
    private static IBISKNIMENodesPlugin IKN_PLUGIN;

    /**
     * The constructor.
     */
    public IBISKNIMENodesPlugin() {
        super();
        IKN_PLUGIN = this;
    }
    
    private Thread t;

    /**
     * This method is called upon plug-in activation.
     * 
     * @param context
     *            The OSGI bundle context
     * @throws Exception
     *             If this GKN_PLUGIN could not be started
     */
    @Override
    public void start(final BundleContext context) throws Exception {
        super.start(context);
        IKN_PLUGIN = this;
       
        IBISKNIMENodesPluginStartupMessageProvider iknpsmp = new IBISKNIMENodesPluginStartupMessageProvider();
        
        String warn_message = "";
        
        //check tool paths
        for(String s: PATHS) {
			String path = getStringPreference(s);
			if(Files.notExists(Paths.get(path)) && !path.equals("")) {
				getDefault().getPreferenceStore().setToDefault(s);
				warn_message += "Path to "+s+" has been changed!"+ System.getProperty("line.separator"); 
				if(s.equals(DB_FILE)) {
					getDefault().getPreferenceStore().setToDefault(USE_HTE);
					getDefault().getPreferenceStore().setToDefault(NOTIFY);
					warn_message += "You can reset the path to hte.db or create a new database file in the KNIME4NGS preference page!" +
							System.getProperty("line.separator"); 
				}
			}
		}
        
        if(warn_message.length()>0) {
        	warn_message += "Please edit the paths in the KNIME4NGS preference page!"+System.getProperty("line.separator");
        }
        
        StartupMessage sm;
        if(warn_message.length()>0) {
            sm = new StartupMessage(warn_message, "Paths changed!",StartupMessage.WARNING, context.getBundle());
        } else {
        	sm = new StartupMessage("KNIME4NGS nodes started correctly.", StartupMessage.INFO, context.getBundle());
        }
        iknpsmp.addMessage(sm);
        
    }
    
	public void startSearchThread(String dir, Table table, Button search, Button cancel) {
		
		t = new Thread(new Runnable() {
		      public void run() {
		    	  
		    	  HashMap<String, String> tool2path = selectSearchDir(dir);
		    	  if(!Thread.currentThread().isInterrupted()) {
		    		  StringBuilder sb = new StringBuilder();
		    		  if(tool2path.size()==0) {
		    			  sb.append("No tools have been found!");
		    		  } else {
		    			  sb.append("The following tools have been found:"+System.getProperty("line.separator"));
		    		  }
		    		  for(String tool: tool2path.keySet()) {
		    			  String path = tool2path.get(tool);
		    			  if(path!=null) {
		    				  sb.append(path+System.getProperty("line.separator"));
		    			  }
		    		  }
		    		  int n = JOptionPane.showConfirmDialog(null,
		    				  sb.toString(),
		    				  "Searching directory "+dir+" finished",
		    				  JOptionPane.OK_CANCEL_OPTION);
		    		  if(n == 0) {
		    			  for(String tool: tool2path.keySet()) {
		    				  String path = tool2path.get(tool);
			    			  if(path!=null) {
			    				  IBISKNIMENodesPlugin.setStringPreference(tool, tool2path.get(tool));
			    			  }
		    			  }
		    		  }
		    		  if(!table.isDisposed()) {
	    				  Display.getDefault().asyncExec(new Runnable () {
							@Override
							public void run() {
								for(TableItem i: table.getItems()) {
									i.setText(1,IBISKNIMENodesPlugin.getStringPreference(i.getText(0)));
								}	
								cancel.setEnabled(false);
								search.setEnabled(true);
							}
	    				  });		    				  
	    			  }
		    	  }
		      }
		    });
		t.start();
	}
	
	public HashMap<String, String> selectSearchDir(String dir){
		if(dir==null) return null;

		HashMap<String, String> tool2path = new HashMap<>();
		
		for(String s: IBISKNIMENodesPlugin.TOOLS.keySet()) {
			if(IBISKNIMENodesPlugin.getStringPreference(s).equals("")) {
				String path = BinaryHandler.checkToolAvailability(s, dir);
				if(path!=null) {
					tool2path.put(s, path);
				}
			}
		}
		return tool2path;
	}
	
	public void cancelSearchThread() {
		if(t.isAlive()) {
			System.out.println("interrupt t");
			t.interrupt();
		}
		System.out.println("thread not alive");
	}
	
	public boolean isSearching() {
		if(t==null) return false;
		return t.isAlive();
	}

    /**
     * This method is called when the plug-in is stopped.
     * 
     * @param context
     *            The OSGI bundle context
     * @throws Exception
     *             If this GKN_PLUGIN could not be stopped
     */
    @Override
    public void stop(final BundleContext context) throws Exception {
    	IKN_PLUGIN = null;
        super.stop(context);
    }

    /**
     * Returns the shared instance.
     * 
     * @return Singleton instance of the Plugin
     */
    public static IBISKNIMENodesPlugin getDefault() {
        return IKN_PLUGIN;
    }

	/**
	 * Initializes a preference store with default preference values 
	 * for this plug-in.
	 * 
	 * @param store the preference store to fill
	 */
	protected void initializeDefaultPreferences(IPreferenceStore store) {
		store.setDefault(USE_HTE, HTE_DEFAULT);
		store.setDefault(THRESHOLD, THRESHOLD_DEFAULT);
		store.setDefault(NOTIFY, NOTIFY_DEFAULT);
		store.setDefault(EMAIL_HOST, EMAIL_HOST_DEFAULT);
		store.setDefault(EMAIL_SENDER, EMAIL_SENDER_DEFAULT);
		store.setDefault(EMAIL_RECEIVER, EMAIL_RECEIVER_DEFAULT);
		
		for(String t: PATHS) {
			store.setDefault(t, "");
		}
	}
	
	public static boolean getBooleanPreference(String id) {
		return IBISKNIMENodesPlugin.getDefault().getPreferenceStore().getBoolean(id);
	}
	
	public static void setBooleanPreference(String id, boolean value) {
		IBISKNIMENodesPlugin.getDefault().getPreferenceStore().setValue(id, value);
	}
	
	public static String getStringPreference(String id) {
		return IBISKNIMENodesPlugin.getDefault().getPreferenceStore().getString(id);
	}
	
	public static void setStringPreference(String id, String value) {
		IBISKNIMENodesPlugin.getDefault().getPreferenceStore().setValue(id, value);
	}
	
	public static void setAllFieldsToDefault() {
		for(String f: FIELDS) {
			IBISKNIMENodesPlugin.getDefault().getPreferenceStore().setToDefault(f);
		}
		for(String s: PATHS) {
			IBISKNIMENodesPlugin.getDefault().getPreferenceStore().setToDefault(s);
		}
	}
}
