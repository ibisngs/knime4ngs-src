package de.helmholtz_muenchen.ibis.knime;

import java.util.ArrayList;
import java.util.List;

import org.knime.workbench.ui.startup.StartupMessage;
import org.knime.workbench.ui.startup.StartupMessageProvider;

public class IBISKNIMENodesPluginStartupMessageProvider implements StartupMessageProvider {


	private static ArrayList<StartupMessage> mylist = new ArrayList<>(); 
	
	public void addMessage(StartupMessage s) {
		mylist.add(s);
	}
	
	@Override
	public List<StartupMessage> getMessages() {
		return mylist;
	}

}
