module.exports = function(grunt) {
	grunt.initConfig({
		uglify: {
			my_target: {
				files: {
					'mesonet.min.js': ['winddata.js', 'windmap.js']
				}
			}
		},
		shell: {
			python: {
				options: {
					stdout: true
				},
				command: 'python interpolation.py'
			}
		}
	});

	grunt.loadNpmTasks('grunt-contrib-uglify');
	grunt.loadNpmTasks('grunt-shell');

	grunt.registerTask('default', ['shell:python', 'uglify']);
};