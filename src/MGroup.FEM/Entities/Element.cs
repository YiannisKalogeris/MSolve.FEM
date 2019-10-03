using System.Collections.Generic;
using System.Linq;
using MGroup.FEM.Interfaces;
using MGroup.MSolve.Discretization.FreedomDegrees;
using MGroup.MSolve.Discretization.Interfaces;

//TODO: Delete this ASAP. 1) Its purpose is element-node connectivity, which should be done through interfaces and inheritance,
//      2) The order of the nodes should be defined by what is now called ElementType
namespace MGroup.FEM.Entities
{
	/// <summary>
	/// Enum defining absorption types.
	/// </summary>
	public enum AbsorptionType
	{
		Unknown = 0,
		Compressional = 1,
		Shear = 2,
	}

	/// <summary>
	/// Class that acts a wrapper element.
	/// </summary>
	public class Element : IElement
	{
		private readonly Dictionary<int, Node> nodesDictionary = new Dictionary<int, Node>();
		private readonly IList<Node> embeddedNodes = new List<Node>();

		/// <summary>
		/// <see cref="Dictionary{TKey,TValue}"/> that contains node ID, node pairs.
		/// </summary>
		public Dictionary<int, Node> NodesDictionary => nodesDictionary;

		/// <summary>
		/// <see cref="Dictionary{TKey,TValue}"/> that contains dof, absorption pairs.
		/// </summary>
		public Dictionary<IDofType, AbsorptionType> Absorptions { get; } = new Dictionary<IDofType, AbsorptionType>();

		/// <summary>
		/// A <see cref="List{T}"/> containing the element nodes.
		/// </summary>
		IReadOnlyList<INode> IElement.Nodes => nodesDictionary.Values.ToList<INode>();

		/// <summary>
		/// A <see cref="List{T}"/> containing the element nodes.
		/// </summary>
		public IList<Node> Nodes => nodesDictionary.Values.ToList();

		/// <summary>
		/// /// <summary>
		/// A <see cref="List{T}"/> containing the nodes embedded in the element.
		/// </summary>
		/// </summary>
		public IList<Node> EmbeddedNodes => embeddedNodes;

		/// <summary>
		/// The element ID.
		/// </summary>
		public int ID { get; set; }

		/// <summary>
		/// The type of the element.
		/// </summary>
		IElementType IElement.ElementType => ElementType;

		/// <summary>
		/// The specific element type.
		/// </summary>
		public IFiniteElement ElementType { get; set; }

		/// <summary>
		/// The <see cref="Subdomain"/> this element belongs to.
		/// </summary>
		ISubdomain IElement.Subdomain => this.Subdomain;

		/// <summary>
		/// The <see cref="Subdomain"/> this element belongs to.
		/// </summary>
		public Subdomain Subdomain { get; set; }

		/// <summary>
		/// Adds a single node to the element.
		/// </summary>
		/// <param name="node">A <see cref="Node"/>.</param>
		public void AddNode(Node node) => nodesDictionary.Add(node.ID, node);

		/// <summary>
		/// Adds a list of nodes to the element.
		/// </summary>
		/// <param name="nodes">A <see cref="List{T}"/> of <see cref="Nodes"/>.</param>
		public void AddNodes(IList<Node> nodes)
		{
			foreach (Node node in nodes) AddNode(node);
		}
	}
}
