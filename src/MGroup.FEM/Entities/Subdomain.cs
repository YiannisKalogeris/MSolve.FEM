using System.Collections.Generic;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.Discretization.Commons;
using MGroup.MSolve.Discretization.FreedomDegrees;
using MGroup.MSolve.Discretization.Interfaces;

namespace MGroup.FEM.Entities
{
	/// <summary>
	/// An entity that defines a subdomain in FEM method.
	/// </summary>
	public class Subdomain : ISubdomain
	{
		private readonly List<Node> nodes = new List<Node>();

		/// <summary>
		/// Defines a <see cref="Subdomain"/>.
		/// </summary>
		/// <param name="id">The Id of the subdomain.</param>
		public Subdomain(int id)
		{
			this.ID = id;
		}

		/// <summary>
		///  A <see cref="Table{TRow,TColumn,TValue}"/> that contains the constrained dofs of the subdomain.
		/// </summary>
		public Table<INode, IDofType, double> Constraints { get; } = new Table<INode, IDofType, double>();

		/// <summary>
		/// A <see cref="List{T}"/> with the elements of the subdomain.
		/// </summary>
		IReadOnlyList<IElement> ISubdomain.Elements => Elements;

		/// <summary>
		/// A <see cref="List{T}"/> with the elements of the subdomain.
		/// </summary>
		public List<Element> Elements { get; } = new List<Element>();

		/// <summary>
		/// The element ID.
		/// </summary>
		public int ID { get; }

		/// <summary>
		/// A <see cref="List{T}"/> with the nodes of the subdomain.
		/// </summary>
		IReadOnlyList<INode> ISubdomain.Nodes => nodes;

		/// <summary>
		/// A <see cref="List{T}"/> with the nodes of the subdomain.
		/// </summary>
		public IReadOnlyList<Node> Nodes => nodes;

		/// <summary>
		/// The ordering of constrained dofs. For more info please refer to <see cref="ISubdomainConstrainedDofOrdering"/>.
		/// </summary>
		public ISubdomainConstrainedDofOrdering ConstrainedDofOrdering { get; set; }

		/// <summary>
		/// The ordering of free dofs. For more info please refer to <see cref="ISubdomainFreeDofOrdering"/>.
		/// </summary>
		public ISubdomainFreeDofOrdering FreeDofOrdering { get; set; }

		/// <summary>
		/// The forces vector of the subdomain.
		/// </summary>
		public Vector Forces { get; set; }

		/// <summary>
		/// Boolean that denotes if the stiffness matrix has been modified.
		/// </summary>
		public bool StiffnessModified { get; set; } = true;

		/// <summary>
		/// Boolean that denotes if the connectivity has been modified.
		/// </summary>
		public bool ConnectivityModified { get; set; } = true;

		/// <summary>
		/// Calculates the elemental displacements.
		/// </summary>
		/// <param name="element">An <see cref="Element"/>.</param>
		/// <param name="constraintScalingFactor">The constraint scaling factor.</param>
		/// <returns>A double aray with the elemental nodal displacements.</returns>
		public double[] CalculateElementIncrementalConstraintDisplacements(IElement element, double constraintScalingFactor)//QUESTION: would it be maybe more clear if we passed the constraintsDictionary as argument??
		{
			var elementNodalDisplacements = new double[FreeDofOrdering.CountElementDofs(element)];
			SubdomainConstrainedDofOrderingBase.ApplyConstraintDisplacements(element, elementNodalDisplacements, Constraints);
			return elementNodalDisplacements;
		}

		/// <summary>
		/// Calculates the elemental displacements.
		/// </summary>
		/// <param name="element">An <see cref="Element"/>.</param>
		/// <param name="globalDisplacementVector">The global displacements vector.</param>
		/// <returns>A double aray with the elemental nodal displacements.</returns>
		public double[] CalculateElementDisplacements(IElement element, IVectorView globalDisplacementVector)//QUESTION: would it be maybe more clear if we passed the constraintsDictionary as argument??
		{
			double[] elementNodalDisplacements = FreeDofOrdering.ExtractVectorElementFromSubdomain(element, globalDisplacementVector);
			SubdomainConstrainedDofOrderingBase.ApplyConstraintDisplacements(element, elementNodalDisplacements, Constraints);
			return elementNodalDisplacements;
		}

		/// <summary>
		/// Clear any saved material stresses of the element.
		/// </summary>
		public void ClearMaterialStresses()
		{
			foreach (var element in Elements) element.ElementType.ClearMaterialStresses();
		}

		/// <summary>
		/// Defines the subdomain nodes from the elements. 
		/// </summary>
		public void DefineNodesFromElements()
		{
			nodes.Clear();
			var nodeComparer = Comparer<Node>.Create((Node node1, Node node2) => node1.ID - node2.ID);
			var nodeSet = new SortedSet<Node>(nodeComparer);
			foreach (Element element in Elements)
			{
				foreach (Node node in element.Nodes) nodeSet.Add(node);
			}

			nodes.AddRange(nodeSet);
		}

		/// <summary>
		/// Extract the subdomain constraints from the global ones.
		/// </summary>
		/// <param name="globalConstraints">A <see cref="Table{TRow,TColumn,TValue}"/> containing the global <see cref="Model"/> constrains.</param>
		public void ExtractConstraintsFromGlobal(Table<INode, IDofType, double> globalConstraints)
		{
			foreach (Node node in Nodes)
			{
				bool isNodeConstrained = globalConstraints.TryGetDataOfRow(node,
					out IReadOnlyDictionary<IDofType, double> constraintsOfNode);
				if (isNodeConstrained)
				{
					foreach (var dofDisplacementPair in constraintsOfNode)
					{
						Constraints[node, dofDisplacementPair.Key] = dofDisplacementPair.Value;
					}
				}
			}
		}

		/// <summary>
		/// Gets the subdomain load vector from the solution.
		/// </summary>
		/// <param name="solution">The solution vector.</param>
		/// <param name="dSolution">The incremental solution vector.</param>
		/// <returns>A vector with the subdomain forces.</returns>
		public IVector GetRhsFromSolution(IVectorView solution, IVectorView dSolution)
		{
			var forces = Vector.CreateZero(FreeDofOrdering.NumFreeDofs); //TODO: use Vector
			foreach (Element element in Elements)
			{
				double[] localSolution = CalculateElementDisplacements(element, solution);
				double[] localdSolution = CalculateElementDisplacements(element, dSolution);
				element.ElementType.CalculateStresses(element, localSolution, localdSolution);
				if (element.ElementType.MaterialModified)
					element.Subdomain.StiffnessModified = true;
				var f = element.ElementType.CalculateForces(element, localSolution, localdSolution);
				FreeDofOrdering.AddVectorElementToSubdomain(element, f, forces);
			}
			return forces;
		}

		/// <summary>
		/// Resets the material modified property of elements.
		/// </summary>
		public void ResetMaterialsModifiedProperty()
		{
			this.StiffnessModified = false;
			foreach (Element element in Elements) element.ElementType.ResetMaterialModified();
		}

		/// <summary>
		/// Saves the materials states of the subdomain elements.
		/// </summary>
		public void SaveMaterialState()
		{
			foreach (Element element in Elements) element.ElementType.SaveMaterialState();
		}

		/// <summary>
		/// Scales the subdomain constraints.
		/// </summary>
		/// <param name="scalingFactor">A double value.</param>
		public void ScaleConstraints(double scalingFactor) => Constraints.ModifyValues((u) => scalingFactor * u);

		/// <summary>
		/// Get the force vector from solution when initial displacements are defined in the model.
		/// </summary>
		/// <param name="solution">Solution vector.</param>
		/// <param name="dSolution">Incremental solution vector.</param>
		/// <param name="boundaryNodes"><see cref="Dictionary{TKey,TValue}"/> containing the boundary nodes.</param>
		/// <param name="initialConvergedBoundaryDisplacements">Initial boundary displacements.</param>
		/// <param name="totalBoundaryDisplacements">Total boundary displacements.</param>
		/// <param name="nIncrement">Increment number.</param>
		/// <param name="totalIncrements">Total number of increments.</param>
		/// <returns>An <see cref="IVector"/> with the subdomain forces.</returns>
		public IVector GetRHSFromSolutionWithInitialDisplacementsEffect(IVectorView solution, IVectorView dSolution, Dictionary<int, INode> boundaryNodes,
			Dictionary<int, Dictionary<IDofType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<IDofType, double>> totalBoundaryDisplacements,
			int nIncrement, int totalIncrements)
		{
			var forces = Vector.CreateZero(FreeDofOrdering.NumFreeDofs); //TODO: use Vector
			foreach (Element element in Elements)
			{
				var localSolution = GetLocalVectorFromGlobalWithoutPrescribedDisplacements(element, solution);
				ImposePrescribedDisplacementsWithInitialConditionSEffect(element, localSolution, boundaryNodes, initialConvergedBoundaryDisplacements, totalBoundaryDisplacements, nIncrement, totalIncrements);
				var localdSolution = GetLocalVectorFromGlobalWithoutPrescribedDisplacements(element, dSolution);
				ImposePrescribed_d_DisplacementsWithInitialConditionSEffect(element, localdSolution, boundaryNodes, initialConvergedBoundaryDisplacements, totalBoundaryDisplacements, nIncrement, totalIncrements);
				element.ElementType.CalculateStresses(element, localSolution, localdSolution);
				if (element.ElementType.MaterialModified)
					element.Subdomain.StiffnessModified = true;
				var f = element.ElementType.CalculateForces(element, localSolution, localdSolution);
				FreeDofOrdering.AddVectorElementToSubdomain(element, f, forces);
			}
			return forces;
		}

		/// <summary>
		/// Gets the elemental nodal displacements from the global model ones.
		/// </summary>
		/// <param name="element">An <see cref="Element"/>.</param>
		/// <param name="globalDisplacementVector">The global displacement vector.</param>
		/// <returns>A double array containing the elemental displacements.</returns>
		public double[] GetLocalVectorFromGlobalWithoutPrescribedDisplacements(Element element, IVectorView globalDisplacementVector)
		{
			double[] elementNodalDisplacements = FreeDofOrdering.ExtractVectorElementFromSubdomain(element, globalDisplacementVector);
			return elementNodalDisplacements;
		}

		/// <summary>
		/// Imposes prescribed displacements.
		/// </summary>
		/// <param name="element">An <see cref="Element"/>.</param>
		/// <param name="localSolution">The solution for the subdomain.</param>
		/// <param name="boundaryNodes">The boundary nodes.</param>
		/// <param name="initialConvergedBoundaryDisplacements">The initial boundary displacements.</param>
		/// <param name="totalBoundaryDisplacements">The total boundary displacements.</param>
		/// <param name="nIncrement">The increment number.</param>
		/// <param name="totalIncrements">The total number of increments.</param>
		public void ImposePrescribedDisplacementsWithInitialConditionSEffect(Element element, double[] localSolution, Dictionary<int, INode> boundaryNodes,
			Dictionary<int, Dictionary<IDofType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<IDofType, double>> totalBoundaryDisplacements,
			int nIncrement, int totalIncrements)
		{

			var elementDOFTypes = element.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(element);
			var matrixAssemblyNodes = element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element);
			int iElementMatrixColumn = 0;
			for (int j = 0; j < elementDOFTypes.Count; j++)
			{
				INode nodeColumn = matrixAssemblyNodes[j];
				int nodalDofsNumber = elementDOFTypes[j].Count;
				if (boundaryNodes.ContainsKey(nodeColumn.ID))
				{
					Dictionary<IDofType, double> nodalConvergedDisplacements = initialConvergedBoundaryDisplacements[nodeColumn.ID];
					Dictionary<IDofType, double> nodalTotalDisplacements = totalBoundaryDisplacements[nodeColumn.ID];
					int positionOfDofInNode = 0;
					foreach (IDofType doftype1 in elementDOFTypes[j])
					{
						if (nodalConvergedDisplacements.ContainsKey(doftype1))
						{
							localSolution[iElementMatrixColumn + positionOfDofInNode] = nodalConvergedDisplacements[doftype1] + (nodalTotalDisplacements[doftype1] - nodalConvergedDisplacements[doftype1]) * ((double)nIncrement / (double)totalIncrements);
							// TODO: this can be done faster: create a dictionary<...,dictionary> with the difference of the two values and use that and precalculate coefficient for scaling
						}
						positionOfDofInNode += 1;
					}
				}
				iElementMatrixColumn += nodalDofsNumber;
			}

		}

		/// <summary>
		/// Imposes prescribed displacements.
		/// </summary>
		/// <param name="element">An <see cref="Element"/>.</param>
		/// <param name="localSolution">The solution for the subdomain.</param>
		/// <param name="boundaryNodes">The boundary nodes.</param>
		/// <param name="initialConvergedBoundaryDisplacements">The initial boundary displacements.</param>
		/// <param name="totalBoundaryDisplacements">The total boundary displacements.</param>
		/// <param name="nIncrement">The increment number.</param>
		/// <param name="totalIncrements">The total number of increments.</param>
		public void ImposePrescribed_d_DisplacementsWithInitialConditionSEffect(Element element, double[] localSolution, Dictionary<int, INode> boundaryNodes,
			Dictionary<int, Dictionary<IDofType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<IDofType, double>> totalBoundaryDisplacements,
			int nIncrement, int totalIncrements)
		{

			var elementDOFTypes = element.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(element);
			var matrixAssemblyNodes = element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element);
			int iElementMatrixColumn = 0;
			for (int j = 0; j < elementDOFTypes.Count; j++)
			{
				INode nodeColumn = matrixAssemblyNodes[j];
				int nodalDofsNumber = elementDOFTypes[j].Count;
				if (boundaryNodes.ContainsKey(nodeColumn.ID))
				{
					Dictionary<IDofType, double> nodalConvergedDisplacements = initialConvergedBoundaryDisplacements[nodeColumn.ID];
					Dictionary<IDofType, double> nodalTotalDisplacements = totalBoundaryDisplacements[nodeColumn.ID];
					int positionOfDofInNode = 0;
					foreach (IDofType doftype1 in elementDOFTypes[j])
					{
						if (nodalConvergedDisplacements.ContainsKey(doftype1))
						{
							localSolution[iElementMatrixColumn + positionOfDofInNode] = (nodalTotalDisplacements[doftype1] - nodalConvergedDisplacements[doftype1]) * ((double)nIncrement / (double)totalIncrements);
							// 1) den vazoume mono (1/increments) alla (nIncrement/increments) dioti metaxu aftwn twn nIncrements den exei mesolavhsei save sta material ths mikroklimakas
							// TODO: this can be done faster: create a dictionary<...,dictionary> with the difference of the two values and use that and precalculate coefficient for scaling
						}
						positionOfDofInNode += 1;
					}
				}
				iElementMatrixColumn += nodalDofsNumber;
			}

		}
	}
}
