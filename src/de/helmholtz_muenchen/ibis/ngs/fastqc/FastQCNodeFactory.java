package de.helmholtz_muenchen.ibis.ngs.fastqc;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "FastQC" Node.
 * 
 *
 * @author Max
 */
public class FastQCNodeFactory 
        extends NodeFactory<FastQCNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public FastQCNodeModel createNodeModel() {
        return new FastQCNodeModel();
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
    public NodeView<FastQCNodeModel> createNodeView(final int viewIndex,
            final FastQCNodeModel nodeModel) {
        return new FastQCNodeView(nodeModel);
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
        return new FastQCNodeDialog();
    }

}

