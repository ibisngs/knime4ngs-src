package de.helmholtz_muenchen.ibis.ngs.bwa;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "BWA" Node.
 * 
 *
 * @author Jan
 */
public class BWANodeFactory 
        extends NodeFactory<BWANodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public BWANodeModel createNodeModel() {
        return new BWANodeModel();
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
    public NodeView<BWANodeModel> createNodeView(final int viewIndex,
            final BWANodeModel nodeModel) {
        return new BWANodeView(nodeModel);
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
        return new BWANodeDialog();
    }

}

