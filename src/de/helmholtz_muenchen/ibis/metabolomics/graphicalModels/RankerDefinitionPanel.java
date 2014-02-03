package de.helmholtz_muenchen.ibis.metabolomics.graphicalModels;

import java.awt.GridLayout;

import javax.swing.BorderFactory;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

public final class RankerDefinitionPanel extends JPanel {

	private static final long serialVersionUID = -7554277935846482740L;
	private final JComboBox<String> m_rankerMethod;
	private final JTextArea m_rankerParam;
	
	
	RankerDefinitionPanel(String name, String ranker, String params){
		super(new GridLayout(0,2));
		this.setBorder(BorderFactory.createTitledBorder(" "+ name + " "));

		// ranker
        JLabel label_rankerMethod = new JLabel("Ranking Method");
        this.add(label_rankerMethod);
        m_rankerMethod = new JComboBox<String>(GraphicalModelsNodeModel.GRAFO_RANKERS);
        m_rankerMethod.setSelectedItem(ranker);
        this.add(m_rankerMethod);
        
        // params
        JLabel label_rankerParams = new JLabel("Ranking Parameters");
        this.add(label_rankerParams);
        m_rankerParam = new JTextArea(10,20);
        m_rankerParam.setText(params);
        JScrollPane rankerParamScrollPane = new JScrollPane(m_rankerParam);
        this.add(rankerParamScrollPane);
	}
	
	public String getRanker(){
		return((String)this.m_rankerMethod.getSelectedItem());
	}
	
	public String getParams(){
		return(m_rankerParam.getText());
	}
}
