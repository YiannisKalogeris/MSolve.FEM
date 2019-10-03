using System;
using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.FEM.Interfaces;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.FreedomDegrees;
using MGroup.MSolve.Discretization.Interfaces;
using MGroup.MSolve.Discretization.Loads;
using MGroup.MSolve.Discretization.Mesh;

namespace MGroup.FEM.Elements
{
	/// <summary>
	/// A concentrated mass.
	/// </summary>
	public class ConcentratedMass3D : IStructuralFiniteElement
	{
		private static readonly IDofType[] nodalDOFTypes = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
		private static readonly IDofType[][] dofs = new IDofType[][] { nodalDOFTypes };
		private readonly double massCoefficient;
		private IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();

		/// <summary>
		/// The element ID.
		/// </summary>
		public int ID => 998;

		/// <summary>
		/// Saves the geometry state of the element.
		/// </summary>
		public CellType CellType { get; } = CellType.Unknown;

		/// <summary>
		/// Dimensionality of the element.
		/// </summary>
		public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;

		/// <summary>
		/// Defines the way that elemental degrees of freedom will be enumerated.
		/// For further info see <see cref="IElementDofEnumerator"/>.
		/// </summary>
		public IElementDofEnumerator DofEnumerator
		{
			get { return dofEnumerator; }
			set { dofEnumerator = value; }
		}

		/// <summary>
		/// Retrieves the dofs of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="ConcentratedMass3D"/>.</param>
		/// <returns>A <see cref="IReadOnlyList{T}"/> that contains a <see cref="IReadOnlyList{T}"/> of <see cref="IDofType"/> with degrees of freedom for each elemental <see cref="Node"/>.</returns>
		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element)
		{
			if (element == null) return dofs;

			var d = new List<List<IDofType>>();
			foreach (var node in element.Nodes)
			{
				var nodeDofs = new List<IDofType>();
				nodeDofs.AddRange(nodalDOFTypes);
				d.Add(nodeDofs);
			}
			return d;
		}

		/// <summary>
		/// Boolean denoting if the material of the element has been modified.
		/// </summary>
		public bool MaterialModified => false;

		/// <summary>
		/// Defines a <see cref="ConcentratedMass3D"/>.
		/// </summary>
		/// <param name="massCoefficient">The mass coefficient.</param>
		public ConcentratedMass3D(double massCoefficient)
		{
			this.massCoefficient = massCoefficient;
		}

		/// <summary>
		/// Defines a <see cref="ConcentratedMass3D"/>.
		/// </summary>
		/// <param name="massCoefficient">The mass coefficient.</param>
		/// <param name="dofEnumerator">An <see cref="IElementDofEnumerator"/>.</param>
		public ConcentratedMass3D(double massCoefficient, IElementDofEnumerator dofEnumerator)
			: this(massCoefficient)
		{
			this.dofEnumerator = dofEnumerator;
		}

		/// <summary>
		/// Calculates the Mass Matrix.
		/// </summary>
		/// <param name="element">An element of type <see cref="ConcentratedMass3D"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the mass matrix of an <see cref="ConcentratedMass3D"/>.</returns>
		public IMatrix MassMatrix(IElement element)
		{
			var mass = Matrix.CreateZero(3, 3);
			mass[0, 0] = massCoefficient;
			mass[1, 1] = massCoefficient;
			mass[2, 2] = massCoefficient;
			return mass;
		}

		/// <summary>
		/// Calculates the stiffness matrix of the element.
		/// </summary>
		/// <param name="element">>An element of type <see cref="ConcentratedMass3D"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the stiffness matrix of an <see cref="ConcentratedMass3D"/>.</returns>
		public IMatrix StiffnessMatrix(IElement element) => Matrix.CreateZero(3, 3);

		/// <summary>
		/// Calculates the damping matrix of the element.
		/// </summary>
		/// <param name="element">>An element of type <see cref="ConcentratedMass3D"/>.</param>
		/// <returns>An <see cref="IMatrix"/> containing the damping matrix of an <see cref="ConcentratedMass3D"/>.</returns>
		public IMatrix DampingMatrix(IElement element) => Matrix.CreateZero(3, 3);

		/// <summary>
		/// Resets any saved material states of the element to its initial state.
		/// </summary>
		public void ResetMaterialModified() { }

		/// <summary>
		/// This method calculates the stresses of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="ConcentratedMass3D"/>.</param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <param name="localdDisplacements">A <see cref="double"/> array containing the displacements change for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="Tuple{T1,T2}"/> of the stresses and strains of the element.</returns>
		public Tuple<double[], double[]> CalculateStresses(IElement element, double[] localDisplacements,
			double[] localdDisplacements)
			=> new Tuple<double[], double[]>(new double[6], new double[6]);

		/// <summary>
		/// This method is used for retrieving the internal forces of the element for logging purposes.
		/// </summary>
		/// <param name="element">An element of type <see cref="ConcentratedMass3D"/>.</param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array containing the forces all degrees of freedom.</returns>
		public double[] CalculateForcesForLogging(IElement element, double[] localDisplacements)
			=> CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);

		/// <summary>
		/// This method calculates the internal forces of the element.
		/// </summary>
		/// <param name="element">An element of type <see cref="ConcentratedMass3D"/>.</param>
		/// <param name="localDisplacements">A <see cref="double"/> array containing the displacements for the degrees of freedom of the element.</param>
		/// <param name="localdDisplacements">A <see cref="double"/> array containing the displacements change for the degrees of freedom of the element.</param>
		/// <returns>A <see cref="double"/> array containing the forces all degrees of freedom.</returns>
		public double[] CalculateForces(IElement element, double[] localDisplacements, double[] localdDisplacements)
			=> new double[6];

		/// <summary>
		/// Calculates the forces applies to an <see cref="ConcentratedMass3D"/> due to <see cref="MassAccelerationLoad"/>.
		/// </summary>
		/// <param name="element">An element of type <see cref="ConcentratedMass3D"/>.</param>
		/// <param name="loads">A list of <see cref="MassAccelerationLoad"/>. For more info see <seealso cref="MassAccelerationLoad"/>.</param>
		/// <returns>A <see cref="double"/> array containing the forces generates due to acceleration for each degree of freedom.</returns>
		public double[] CalculateAccelerationForces(IElement element, IList<MassAccelerationLoad> loads)
		{
			var accelerations = new double[3];
			IMatrix massMatrix = MassMatrix(element);

			foreach (MassAccelerationLoad load in loads)
			{
				int index = 0;
				foreach (IDofType[] nodalDOFTypes in dofs)
				{
					foreach (IDofType dofType in nodalDOFTypes)
					{
						if (dofType == load.DOF) accelerations[index] += load.Amount;
						index++;
					}
				}
			}

			return massMatrix.Multiply(accelerations);
		}

		/// <summary>
		/// Clear the material state of the element.
		/// </summary>
		public void ClearMaterialState() { }

		/// <summary>
		/// Save the current material state of the element.
		/// </summary>
		public void SaveMaterialState() { }

		/// <summary>
		/// Clear any saved material stresses of the element.
		/// </summary>
		public void ClearMaterialStresses() { }
	}
}
