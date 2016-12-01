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
package de.helmholtz_muenchen.ibis.ngs.pindel;


import java.io.File;
import java.nio.file.Paths;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.node.ExecutionContext;

import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.lofs.FileInputReader;
import de.helmholtz_muenchen.ibis.utils.lofs.PathProcessor;
import de.helmholtz_muenchen.ibis.utils.lofs.Writer_Output;

/**
 * 
 * @author Marie-Sophie Friedl
 * @author Maximilian Hastreiter
 *
 */

public class RunPindel extends PindelNodeModel {
	
	public RunPindel(){
		
	}
	
	
	public String createOutputFilePindel(String filebase, String toolextension, String variant) throws Exception{
			
		return filebase+"."+toolextension;
	}
	
	public void PindelConfig(String bamfile, String ismfile, String configfile) throws Exception {
		
		PindelNodeModel.logger.info("Creating Pindel config file: "+configfile);
		
		// read insert size metrics file -> extract mean insert size
		FileInputReader fr = new FileInputReader(ismfile);		
		String line;
		boolean loopbreak=false;
		while((line=fr.read())!=null){
			if(line.startsWith("MEDIAN_INSERT_SIZE")){
				loopbreak=true;
				break;
			}
		}
		String mean_is="";
		if(loopbreak){
			line=fr.read();			
			String split[] = line.split("\t");
			if(split.length<5){
				throw new Exception("Format error: Insert size metrics file does not contain mean insert size");
			}
			mean_is=split[4];
		}
		else{
			throw new Exception("Format error: File does not contain insert size statistics");
		}
		fr.closer();
		if(mean_is.equals("")){
			throw new Exception("Format error:  Cannot find mean insert size");
		}
		PindelNodeModel.logger.info("Mean insert size: "+mean_is);
		
		/**
		 * Convert double to int to avoid wrong IDs in VCF File
		 */
		int mean = (int)Double.parseDouble(mean_is);
		
		String samplename=Paths.get(PathProcessor.getBase(bamfile)).getFileName().toString();
		
		/* format of config file
		 * bamfile \t mean insert size \t sample name
		 */
		
		// write config file
		Writer_Output wo = new Writer_Output(configfile);
		wo.writeFile(bamfile);
		wo.writeFile("\t");
		wo.writeFile(mean+"");
		wo.writeFile("\t");
		wo.writeFile(samplename);	
		
		wo.closew();
		
		PindelNodeModel.logger.info("Config file can be found in "+configfile);
		
	}
	
	
	public void Pindel (ExecutionContext exec, String pex, String conf, String ref, String out, String intv, int [] threadsbins, double [] params) throws Exception{
		
		String cmd = pex;
		cmd+=" -i "+conf;
		cmd+=" -f "+ref;
		cmd+=" -o "+out;
		
		if(!intv.equals("")){
			cmd+=" -c "+intv;
		}
		else{
			cmd+=" -c ALL";
		}
		
		cmd+=" -T "+threadsbins[0];
		cmd+=" -w "+threadsbins[1];
		
		cmd+=" -d "+(int) params[0];
		cmd+=" -a "+(int) params[1];
		cmd+=" -m "+(int) params[2];
		cmd+=" -u "+params[3];
		cmd+=" -e "+params[4];
		
		PindelNodeModel.logger.info("Running Pindel...");
		PindelNodeModel.logger.info("Log files can be found in "+out+".out.log and "+out+".err.log");
		
    	/**Execute**/
    	String lockFile = out + SuccessfulRunChecker.LOCK_ENDING;
    	super.executeCommand(new String[]{StringUtils.join(cmd, " ")}, out, exec, new File(lockFile),out+".out.log", out+".err.log");
			
	}
	
	public void Pindel2VCF (ExecutionContext exec, String p2vcf, String ref, String refn, String refd, String pfile, String vout, boolean [] fs, double [] nums) throws Exception{
		
		String cmd = p2vcf;
		cmd+=" -r "+ref;
		cmd+=" -R "+refn;
		cmd+=" -d "+refd;
		cmd+=" -p "+pfile;
		cmd+=" -v "+vout;
		
		if(fs[0]){
			cmd+=" -G";
		}
		if(fs[1]){
			cmd+=" -b";
		}
		
		cmd+=" -mc "+(int)nums[0];
		cmd+=" -he "+nums[1];
		cmd+=" -ho "+nums[2];
		cmd+=" -e "+nums[3];
		cmd+=" -is "+nums[4];
		
		if(fs[2]){
			cmd+=" -as "+nums[5];
		}
		
		PindelNodeModel.logger.info("Running Pindel2VCF...");
		PindelNodeModel.logger.info("Log files can be found in "+vout+".out.log and "+vout+".err.log");
		
    	/**Execute**/
    	String lockFile = vout + SuccessfulRunChecker.LOCK_ENDING;
    	super.executeCommand(new String[]{StringUtils.join(cmd, " ")}, vout, exec, new File(lockFile),vout+".out.log", vout+".err.log");
		
	}

}
