package de.helmholtz_muenchen.ibis.ngs.caseControlAnalyzer;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "CaseControlAnalyzer" Node.
 * 
 *
 * @author Tim Jeske
 */
public class CaseControlAnalyzerNodeFactory 
        extends NodeFactory<CaseControlAnalyzerNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public CaseControlAnalyzerNodeModel createNodeModel() {
        return new CaseControlAnalyzerNodeModel();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getNrNodeViews() {
        return 1;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeView<CaseControlAnalyzerNodeModel> createNodeView(final int viewIndex,
            final CaseControlAnalyzerNodeModel nodeModel) {
        return new CaseControlAnalyzerNodeView(nodeModel);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean hasDialog() {
        return true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeDialogPane createNodeDialogPane() {
        return new CaseControlAnalyzerNodeDialog();
    }

}

