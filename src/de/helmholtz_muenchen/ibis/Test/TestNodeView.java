package de.helmholtz_muenchen.ibis.Test;

import org.knime.core.node.NodeView;

/**
 * <code>NodeView</code> for the "Test" Node.
 * Test node for parallel execution / cluster problem.
 *
 * @author Michael Kluge
 */
public class TestNodeView extends NodeView<TestNodeModel> {

    /**
     * Creates a new view.
     * 
     * @param nodeModel The model (class: {@link TestNodeModel})
     */
    protected TestNodeView(final TestNodeModel nodeModel) {
        super(nodeModel);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void modelChanged() {

        TestNodeModel nodeModel = 
            (TestNodeModel)getNodeModel();
        assert nodeModel != null;
        
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void onClose() {

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void onOpen() {

    }

}

