package de.helmholtz_muenchen.ibis.ngs.summarizelof;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "SummarizeLOF" Node.
 * 
 *
 * @author 
 */
public class SummarizeLOFNodeFactory 
        extends NodeFactory<SummarizeLOFNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public SummarizeLOFNodeModel createNodeModel() {
        return new SummarizeLOFNodeModel();
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
    public NodeView<SummarizeLOFNodeModel> createNodeView(final int viewIndex,
            final SummarizeLOFNodeModel nodeModel) {
        return new SummarizeLOFNodeView(nodeModel);
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
        return new SummarizeLOFNodeDialog();
    }

}

