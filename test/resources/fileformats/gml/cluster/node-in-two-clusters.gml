graph [
	node [
		id 1
	]
	node [
		id 2
	]
	node [
		id 3
	]
	edge [
		source 1
		target 2
	]
	edge [
		source 2
		target 3
	]
	edge [
		source 3
		target 1
	]
]
rootcluster [
	cluster [
		id	1
		vertex "2"
		vertex "3"
	]
	vertex "2"
]
