using System;
using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.MSolve.Discretization.Interfaces;

namespace MGroup.FEM.Interfaces
{
	/// <summary>
	/// An interface for all specific element types.
	/// </summary>
	public interface IFiniteElement : IElementType
	{
		/// <summary>
		/// The element ID.
		/// </summary>
		int ID { get; }

		/// <summary>
		/// The dimensionality of the element.
		/// </summary>
		ElementDimensions ElementDimensions { get; }
	}
}
