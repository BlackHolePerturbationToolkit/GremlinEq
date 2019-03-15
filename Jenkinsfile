pipeline {
  agent { docker { image 'ubuntu-doxygen' } }
  stages {
    stage('Generate documentation') {
      steps {
        withCredentials([sshUserPrivateKey(credentialsId: 'gremlin', keyFileVariable: 'SSH_KEY')]) {
          withEnv(["GIT_SSH_COMMAND=ssh -o StrictHostKeyChecking=no -i ${SSH_KEY}"]) {
            sh 'echo ${SSH_KEY}'
            sh 'echo ${GIT_SSH_COMMAND}'
            dir(path: 'gh-pages') {
              git url: 'git@github.com:BlackHolePerturbationToolkit/GremlinEq.git', branch: 'gh-pages'
            }
            sh '( cat Doxyfile; echo "OUTPUT_DIRECTORY=gh-pages"; echo "HTML_OUTPUT=doc"; echo "GENERATE_LATEX=NO"; echo "SHORT_NAMES=YES";) | doxygen -'
            dir(path: 'gh-pages') {
              sh 'git commit -m "Update documentation"'
              sh 'git push'
            }
          }
        }
      }
    }
  }
}