#extension GL_EXT_gpu_shader4 : enable

varying vec3 fragment_normal;
uniform vec3 L0_pos;
uniform vec3 L1_pos;

float LightIntensity()
{
	//vec3 L0_pos = vec3(1.0, 0.0, 0.0);
	//vec3 L1_pos = vec3(0.0, 0.0, 1.0);
	float L0 = 0.5*max(dot(fragment_normal, L0_pos), 0.0);
	float L1 = 0.5*max(dot(fragment_normal, L1_pos), 0.0);
	
	return L0 + L1;
	//return 0.5;

}

void main(void)
{
	//float error = (1.0 - gl_TexCoord[0].x)*(240.0/360.0);
	float light = LightIntensity();
	vec4 color = gl_Color;//vec4(1.0, 1.0, 1.0, 1.0);
	gl_FragColor = color*light;
	
	//assume that the error is the H component of an HSV triple
	//vec4 HSV = vec4((240.0/60.0)*(1.0 - error), 1.0, 0.5, 1.0);
	//vec4 RGB = HSVtoRGB(HSV);
	
	//convert from HSL to RGB
	//vec4 HSL = vec4(error, 1.0, light, 1.0);
	//vec4 RGB = HSLtoRGB(HSL);
	
	//gl_FragColor = vec4(0.5, 0.5, 0.5, 1.0);
}