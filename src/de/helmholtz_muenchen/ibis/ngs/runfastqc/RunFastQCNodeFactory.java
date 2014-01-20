package de.helmholtz_muenchen.ibis.ngs.runfastqc;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "RunFastQC" Node.
 * 
 *
 * @author 
 */
public class RunFastQCNodeFactory 
        extends NodeFactory<RunFastQCNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public RunFastQCNodeModel createNodeModel() {
        return new RunFastQCNodeModel();
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
    public NodeView<RunFastQCNodeModel> createNodeView(final int viewIndex,
            final RunFastQCNodeModel nodeModel) {
        return new RunFastQCNodeView(nodeModel);
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
        return new RunFastQCNodeDialog();
    }

}

