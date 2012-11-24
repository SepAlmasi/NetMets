varying vec3 fragment_normal;

void main(void)
{
	//fragment_normal = normalize(gl_NormalMatrix * gl_Normal);
	//fragment_normal = vec3(0.0, 0.0, 0.5);
	fragment_normal = gl_Normal;
	gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
	gl_TexCoord[0] = gl_MultiTexCoord0;
	gl_FrontColor = gl_Color;
}
