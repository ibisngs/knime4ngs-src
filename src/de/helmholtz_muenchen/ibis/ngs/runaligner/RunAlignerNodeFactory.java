package de.helmholtz_muenchen.ibis.ngs.runaligner;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "RunAligner" Node.
 * 
 *
 * @author 
 */
public class RunAlignerNodeFactory 
        extends NodeFactory<RunAlignerNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public RunAlignerNodeModel createNodeModel() {
        return new RunAlignerNodeModel();
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
    public NodeView<RunAlignerNodeModel> createNodeView(final int viewIndex,
            final RunAlignerNodeModel nodeModel) {
        return new RunAlignerNodeView(nodeModel);
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
        return new RunAlignerNodeDialog();
    }

}

