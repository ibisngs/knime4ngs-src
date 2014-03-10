package de.helmholtz_muenchen.ibis.Test;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "Test" Node.
 * Test node for parallel execution / cluster problem.
 *
 * @author Michael Kluge
 */
public class TestNodeFactory 
        extends NodeFactory<TestNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public TestNodeModel createNodeModel() {
        return new TestNodeModel();
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
    public NodeView<TestNodeModel> createNodeView(final int viewIndex,
            final TestNodeModel nodeModel) {
        return new TestNodeView(nodeModel);
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
        return new TestNodeDialog();
    }

}

